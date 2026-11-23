use crate::binread::Endian;
use crate::error::{NmrIoError, Result};
use crate::virtual_fs::LoadedDataset;
use nmr_core::{AcquisitionParams, AxisMeta, Domain, NmrFormat, Signal1D, Signal2D};
use num_complex::Complex64;
use std::collections::BTreeMap;

const FDATA_FLOATS: usize = 512;
const FDATA_BYTES: usize = FDATA_FLOATS * 4;

// Key NMRPipe fdata.h indices (verified against fdatap.h)
const IDX_FDMAGIC: usize = 0;
const IDX_FDFLTFORMAT: usize = 1;
const IDX_FDFLTORDER: usize = 2;
const IDX_FDDIMCOUNT: usize = 9;
const IDX_FDF2LABEL: usize = 16;
const IDX_FDF1LABEL: usize = 18;
const IDX_FDF2QUADFLAG: usize = 56;
const IDX_FDF1QUADFLAG: usize = 55;
const IDX_FDF2CAR: usize = 66;
const IDX_FDF1CAR: usize = 67;
const IDX_FDSIZE: usize = 99;
const IDX_FDF2SW: usize = 100;
const IDX_FDF2ORIG: usize = 101;
const IDX_FDQUADFLAG: usize = 106;
const IDX_FDF2OBS: usize = 119;
const IDX_FDF1OBS: usize = 218;
const IDX_FDSPECNUM: usize = 219;
const IDX_FDF2FTFLAG: usize = 220;
const IDX_FDTRANSPOSED: usize = 221;
const IDX_FDF1FTFLAG: usize = 222;
const IDX_FDF1SW: usize = 229;
const IDX_FDF1ORIG: usize = 249;
const IDX_FDTITLE: usize = 297;

/// IEEE little-endian float marker written by NMRPipe into `FDFLTFORMAT`.
const PIPE_FLTFORMAT_BYTES: [u8; 4] = [0xef, 0xee, 0x6e, 0x4f];

/// Quick heuristic check to tell whether `bytes` starts with a 2048-byte NMRPipe `fdata` header.
pub fn looks_like_nmrpipe(bytes: &[u8]) -> bool {
    if bytes.len() < FDATA_BYTES {
        return false;
    }
    if let Some((endian, fdata)) = detect_pipe_endian(bytes) {
        let _ = endian;
        let ndim = fdata[IDX_FDDIMCOUNT];
        let size = fdata[IDX_FDSIZE];
        ndim >= 1.0 && ndim <= 4.0 && size >= 1.0 && size.fract() == 0.0
    } else {
        false
    }
}

fn detect_pipe_endian(bytes: &[u8]) -> Option<(Endian, [f32; FDATA_FLOATS])> {
    if bytes.len() < FDATA_BYTES {
        return None;
    }
    for endian in [Endian::Little, Endian::Big] {
        let mut fdata = [0f32; FDATA_FLOATS];
        for (i, slot) in fdata.iter_mut().enumerate() {
            let b: [u8; 4] = bytes[i * 4..(i + 1) * 4].try_into().unwrap();
            *slot = match endian {
                Endian::Little => f32::from_le_bytes(b),
                Endian::Big => f32::from_be_bytes(b),
            };
        }
        // Standard NMRPipe magic check: FDMAGIC == 0.0 and FDFLTFORMAT matches magic or FDFLTORDER ~ 2.345
        let magic_ok = fdata[IDX_FDMAGIC] == 0.0
            && (fdata[IDX_FDFLTFORMAT].to_bits() == u32::from_le_bytes(PIPE_FLTFORMAT_BYTES)
                || fdata[IDX_FDFLTFORMAT] == 0.0
                || (fdata[IDX_FDFLTORDER] - 2.345).abs() < 0.01);
        let ndim = fdata[IDX_FDDIMCOUNT];
        let size = fdata[IDX_FDSIZE];
        if magic_ok && (1.0..=4.0).contains(&ndim) && ndim.fract() == 0.0 && size >= 1.0 {
            return Some((endian, fdata));
        }
    }
    None
}

fn read_fdata_str(bytes: &[u8], start_idx: usize, count_floats: usize) -> String {
    let start = start_idx * 4;
    let end = (start + count_floats * 4).min(bytes.len());
    let slice = &bytes[start..end];
    let nul = slice.iter().position(|&b| b == 0).unwrap_or(slice.len());
    String::from_utf8_lossy(&slice[..nul]).trim().to_string()
}

/// Load a 1D or 2D NMRPipe binary buffer (`.fid`, `.ft1`, `.ft2`).
pub fn load_nmrpipe_bytes(filename: &str, bytes: &[u8]) -> Result<LoadedDataset> {
    let (endian, fdata) = detect_pipe_endian(bytes).ok_or_else(|| NmrIoError::Malformed {
        format: "NMRPipe",
        detail: "Invalid or unsupported 2048-byte fdata header".to_string(),
    })?;

    let ndim = fdata[IDX_FDDIMCOUNT].round() as usize;
    let size = (fdata[IDX_FDSIZE].round() as usize).max(1);
    let f2_quad_real_only = fdata[IDX_FDF2QUADFLAG] == 1.0;
    let f1_quad_real_only = fdata[IDX_FDF1QUADFLAG] == 1.0;
    let transposed = fdata[IDX_FDTRANSPOSED] == 1.0;
    let f2_ft = fdata[IDX_FDF2FTFLAG] != 0.0;
    let f1_ft = fdata[IDX_FDF1FTFLAG] != 0.0;

    let sf2 = if fdata[IDX_FDF2OBS] > 0.0 {
        fdata[IDX_FDF2OBS] as f64
    } else {
        500.0
    };
    let sw2 = if fdata[IDX_FDF2SW] > 0.0 {
        fdata[IDX_FDF2SW] as f64
    } else {
        6000.0
    };
    let car2 = fdata[IDX_FDF2CAR] as f64;
    let orig2 = fdata[IDX_FDF2ORIG] as f64;
    let offset2_ppm = if sw2 > 0.0 && sf2 > 0.0 && orig2 != 0.0 {
        Some((orig2 + sw2) / sf2)
    } else {
        None
    };

    let f2_label = read_fdata_str(bytes, IDX_FDF2LABEL, 2);
    let f1_label = read_fdata_str(bytes, IDX_FDF1LABEL, 2);
    let title = read_fdata_str(bytes, IDX_FDTITLE, 15);

    let mut vendor_params = BTreeMap::new();
    vendor_params.insert("FDDIMCOUNT".into(), ndim.to_string());
    vendor_params.insert("FDSIZE".into(), size.to_string());
    vendor_params.insert("FDF2SW".into(), sw2.to_string());
    vendor_params.insert("FDF2OBS".into(), sf2.to_string());

    let indirect_axis = if ndim >= 2 {
        let sf1 = if fdata[IDX_FDF1OBS] > 0.0 {
            fdata[IDX_FDF1OBS] as f64
        } else {
            sf2
        };
        let sw1 = if fdata[IDX_FDF1SW] > 0.0 {
            fdata[IDX_FDF1SW] as f64
        } else {
            sw2
        };
        let orig1 = fdata[IDX_FDF1ORIG] as f64;
        Some(AxisMeta {
            nucleus: if f1_label.is_empty() { "1H".into() } else { f1_label },
            spectrometer_freq_mhz: sf1,
            spectral_width_hz: sw1,
            offset_ppm: if orig1 != 0.0 { Some((orig1 + sw1) / sf1) } else { None },
            carrier_ppm: fdata[IDX_FDF1CAR] as f64,
        })
    } else {
        None
    };

    let meta = AcquisitionParams {
        title: if title.is_empty() {
            filename.to_string()
        } else {
            title
        },
        format: NmrFormat::NmrPipe,
        axis: AxisMeta {
            nucleus: if f2_label.is_empty() { "1H".into() } else { f2_label },
            spectrometer_freq_mhz: sf2,
            spectral_width_hz: sw2,
            offset_ppm: offset2_ppm,
            carrier_ppm: if car2 != 0.0 { car2 } else { 5.0 },
        },
        indirect_axis,
        group_delay: None,
        num_scans: 1,
        vendor_params,
    };

    // Read all remaining float32 values after 2048-byte header
    let payload = &bytes[FDATA_BYTES..];
    let n_floats = payload.len() / 4;
    let mut vals = Vec::with_capacity(n_floats);
    for i in 0..n_floats {
        let b: [u8; 4] = payload[i * 4..(i + 1) * 4].try_into().unwrap();
        let v = match endian {
            Endian::Little => f32::from_le_bytes(b),
            Endian::Big => f32::from_be_bytes(b),
        };
        vals.push(v as f64);
    }

    if ndim == 1 {
        let domain = if f2_ft { Domain::Frequency } else { Domain::Time };
        let data: Vec<Complex64> = if f2_quad_real_only {
            vals.into_iter().take(size).map(|r| Complex64::new(r, 0.0)).collect()
        } else if vals.len() >= size * 2 {
            // NMRPipe unappended format: [Re_0..Re_{N-1}, Im_0..Im_{N-1}]
            (0..size)
                .map(|i| Complex64::new(vals[i], vals[size + i]))
                .collect()
        } else {
            vals.into_iter().map(|r| Complex64::new(r, 0.0)).collect()
        };
        Ok(LoadedDataset::D1(Signal1D::new(data, domain, meta)))
    } else {
        let is_real_only = (transposed && f1_quad_real_only) || (!transposed && f2_quad_real_only);
        let multi = if is_real_only { 1 } else { 2 };
        let mut rows = (fdata[IDX_FDSPECNUM].round() as usize).max(1);
        if fdata[IDX_FDQUADFLAG] == 0.0 && is_real_only {
            rows *= 2;
        }
        let stride = size * multi;
        let actual_rows = rows.min(vals.len() / stride.max(1)).max(1);
        let mut data = Vec::with_capacity(actual_rows * size);
        for r in 0..actual_rows {
            let row_slice = &vals[r * stride..((r + 1) * stride).min(vals.len())];
            if is_real_only {
                for &re in row_slice.iter().take(size) {
                    data.push(Complex64::new(re, 0.0));
                }
            } else {
                for c in 0..size {
                    let re = row_slice.get(c).copied().unwrap_or(0.0);
                    let im = row_slice.get(size + c).copied().unwrap_or(0.0);
                    data.push(Complex64::new(re, im));
                }
            }
        }
        Ok(LoadedDataset::D2(Signal2D::new(
            actual_rows,
            size,
            data,
            if f2_ft { Domain::Frequency } else { Domain::Time },
            if f1_ft { Domain::Frequency } else { Domain::Time },
            meta,
        )))
    }
}

/// Serialize a [`Signal1D`] into standard little-endian NMRPipe format bytes.
pub fn write_nmrpipe_1d(sig: &Signal1D) -> Vec<u8> {
    let mut fdata = [0f32; FDATA_FLOATS];
    fdata[IDX_FDMAGIC] = 0.0;
    fdata[IDX_FDFLTFORMAT] = f32::from_le_bytes(PIPE_FLTFORMAT_BYTES);
    fdata[IDX_FDFLTORDER] = 2.345;
    fdata[IDX_FDDIMCOUNT] = 1.0;
    fdata[IDX_FDSIZE] = sig.len() as f32;
    fdata[IDX_FDSPECNUM] = 1.0;
    fdata[IDX_FDQUADFLAG] = 0.0;
    fdata[IDX_FDF2QUADFLAG] = 0.0; // complex
    fdata[IDX_FDF2FTFLAG] = if sig.domain == Domain::Frequency { 1.0 } else { 0.0 };
    fdata[IDX_FDF2SW] = sig.meta.axis.spectral_width_hz as f32;
    fdata[IDX_FDF2OBS] = sig.meta.axis.spectrometer_freq_mhz as f32;
    fdata[IDX_FDF2CAR] = sig.meta.axis.carrier_ppm as f32;
    fdata[IDX_FDF2ORIG] = ((sig.meta.axis.left_ppm() * sig.meta.axis.spectrometer_freq_mhz)
        - sig.meta.axis.spectral_width_hz) as f32;

    let mut out = Vec::with_capacity(FDATA_BYTES + sig.len() * 8);
    for f in fdata {
        out.extend_from_slice(&f.to_le_bytes());
    }
    // Write string fields directly into header bytes
    write_hdr_str(&mut out, IDX_FDF2LABEL, 8, &sig.meta.axis.nucleus);
    write_hdr_str(&mut out, IDX_FDTITLE, 60, &sig.meta.title);

    // Unappended complex layout: Re[0..N], then Im[0..N]
    for c in &sig.data {
        out.extend_from_slice(&(c.re as f32).to_le_bytes());
    }
    for c in &sig.data {
        out.extend_from_slice(&(c.im as f32).to_le_bytes());
    }
    out
}

fn write_hdr_str(buf: &mut [u8], idx_float: usize, max_bytes: usize, text: &str) {
    let start = idx_float * 4;
    let raw = text.as_bytes();
    let n = raw.len().min(max_bytes);
    buf[start..start + n].copy_from_slice(&raw[..n]);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_nmrpipe_1d() {
        let mut meta = AcquisitionParams::default();
        meta.title = "PipeTest".into();
        let sig = Signal1D::new(
            vec![
                Complex64::new(1.0, -2.0),
                Complex64::new(3.5, 4.5),
                Complex64::new(-0.25, 0.75),
            ],
            Domain::Frequency,
            meta,
        );
        let bytes = write_nmrpipe_1d(&sig);
        assert!(looks_like_nmrpipe(&bytes));
        let loaded = load_nmrpipe_bytes("test.ft1", &bytes).unwrap().into_1d().unwrap();
        assert_eq!(loaded.len(), 3);
        assert_eq!(loaded.domain, Domain::Frequency);
        assert_eq!(loaded.meta.title, "PipeTest");
        assert!((loaded.data[1].re - 3.5).abs() < 1e-6);
        assert!((loaded.data[1].im - 4.5).abs() < 1e-6);
    }
}
