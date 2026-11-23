use crate::binread::{BinReader, Endian};
use crate::error::{NmrIoError, Result};
use nmr_core::{AcquisitionParams, AxisMeta, Domain, NmrFormat, Signal2D};
use num_complex::Complex64;
use std::collections::BTreeMap;

/// Parse a 2D Sparky UCSF (`.ucsf`) tiled binary file in pure Rust.
///
/// Big-Endian layout: 180-byte file header + `naxis` × 128-byte axis headers + tiled `f32` matrix.
pub fn load_ucsf_bytes(filename: &str, bytes: &[u8]) -> Result<Signal2D> {
    if bytes.len() < 180 {
        return Err(NmrIoError::Malformed {
            format: "Sparky UCSF",
            detail: "File shorter than 180-byte header".to_string(),
        });
    }

    let mut r = BinReader::new(bytes, Endian::Big, "Sparky UCSF");
    let ident = r.read_ascii_trimmed(10)?;
    if !ident.starts_with("UCSF NMR") {
        return Err(NmrIoError::Malformed {
            format: "Sparky UCSF",
            detail: format!("Unexpected ident '{ident}'"),
        });
    }
    let naxis = r.read_u8()? as usize;
    let _ncomponents = r.read_u8()?;
    let _encoding = r.read_u8()?;
    let _version = r.read_u8()?;
    let _owner = r.read_ascii_trimmed(9)?;
    let _date = r.read_ascii_trimmed(26)?;
    let comment = r.read_ascii_trimmed(80)?;
    r.skip(3)?; // alignment pad before seek_pos
    let _seek_pos = r.read_i32()?;
    r.skip(44)?; // scratch(40) + pad(4)

    if naxis != 2 {
        return Err(NmrIoError::UnsupportedFormat(format!(
            "Sparky UCSF currently expects 2D (found {naxis} axes)"
        )));
    }

    struct UcsfAxis {
        nucleus: String,
        npoints: usize,
        bsize: usize,
        sf: f64,
        sw: f64,
        xmtr_ppm: f64,
    }

    let mut axes = Vec::with_capacity(naxis);
    for _ in 0..naxis {
        let nucleus = r.read_ascii_trimmed(6)?;
        let _spectral_shift = r.read_i16()?;
        let npoints = r.read_u32()? as usize;
        let _size = r.read_u32()? as usize;
        let bsize = (r.read_u32()? as usize).max(1);
        let sf = r.read_f32()? as f64;
        let sw = r.read_f32()? as f64;
        let xmtr = r.read_f32()? as f64;
        let _p0 = r.read_f32()?;
        let _p1 = r.read_f32()?;
        let _first_pt_scale = r.read_f32()?;
        r.skip(84)?;
        axes.push(UcsfAxis {
            nucleus,
            npoints,
            bsize,
            sf: if sf > 0.0 { sf } else { 500.0 },
            sw: if sw > 0.0 { sw } else { 6000.0 },
            xmtr_ppm: xmtr,
        });
    }

    let rows = axes[0].npoints.max(1);
    let cols = axes[1].npoints.max(1);
    let brow = axes[0].bsize;
    let bcol = axes[1].bsize;
    let trows = rows.div_ceil(brow);
    let tcols = cols.div_ceil(bcol);

    let mut matrix = vec![Complex64::new(0.0, 0.0); rows * cols];
    for tr in 0..trows {
        for tc in 0..tcols {
            for r_in_tile in 0..brow {
                for c_in_tile in 0..bcol {
                    if r.pos() + 4 > bytes.len() {
                        break;
                    }
                    let val = r.read_f32()? as f64;
                    let gr = tr * brow + r_in_tile;
                    let gc = tc * bcol + c_in_tile;
                    if gr < rows && gc < cols {
                        matrix[gr * cols + gc] = Complex64::new(val, 0.0);
                    }
                }
            }
        }
    }

    let meta = AcquisitionParams {
        title: if comment.is_empty() {
            filename.to_string()
        } else {
            comment
        },
        format: NmrFormat::SparkyUcsf,
        axis: AxisMeta {
            nucleus: axes[1].nucleus.clone(),
            spectrometer_freq_mhz: axes[1].sf,
            spectral_width_hz: axes[1].sw,
            offset_ppm: None,
            carrier_ppm: axes[1].xmtr_ppm,
        },
        indirect_axis: Some(AxisMeta {
            nucleus: axes[0].nucleus.clone(),
            spectrometer_freq_mhz: axes[0].sf,
            spectral_width_hz: axes[0].sw,
            offset_ppm: None,
            carrier_ppm: axes[0].xmtr_ppm,
        }),
        group_delay: None,
        num_scans: 1,
        vendor_params: BTreeMap::new(),
    };

    Ok(Signal2D::new(
        rows,
        cols,
        matrix,
        Domain::Frequency,
        Domain::Frequency,
        meta,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loads_sparky_ucsf_2d() {
        let mut buf = vec![0u8; 180 + 2 * 128 + 4 * 4];
        buf[0..8].copy_from_slice(b"UCSF NMR");
        buf[10] = 2; // naxis
        buf[11] = 1; // ncomponents
        buf[12] = 0; // encoding
        buf[13] = 2; // version
        // Axis 0 at 180: 2x2 matrix with 2x2 tile
        buf[180..182].copy_from_slice(b"1H");
        buf[188..192].copy_from_slice(&2u32.to_be_bytes());
        buf[192..196].copy_from_slice(&2u32.to_be_bytes());
        buf[196..200].copy_from_slice(&2u32.to_be_bytes());
        buf[200..204].copy_from_slice(&500.0f32.to_be_bytes());
        buf[204..208].copy_from_slice(&5000.0f32.to_be_bytes());
        // Axis 1 at 308
        buf[308..310].copy_from_slice(b"1H");
        buf[316..320].copy_from_slice(&2u32.to_be_bytes());
        buf[320..324].copy_from_slice(&2u32.to_be_bytes());
        buf[324..328].copy_from_slice(&2u32.to_be_bytes());
        buf[328..332].copy_from_slice(&500.0f32.to_be_bytes());
        buf[332..336].copy_from_slice(&5000.0f32.to_be_bytes());
        // 4 floats at 436
        for (i, v) in [1.0f32, 2.0, 3.0, 4.0].iter().enumerate() {
            let p = 436 + i * 4;
            buf[p..p + 4].copy_from_slice(&v.to_be_bytes());
        }
        let s2 = load_ucsf_bytes("cosy.ucsf", &buf).unwrap();
        assert_eq!(s2.rows, 2);
        assert_eq!(s2.cols, 2);
        assert!((s2.data[3].re - 4.0).abs() < 1e-6);
    }
}
