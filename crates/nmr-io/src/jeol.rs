use crate::binread::{BinReader, Endian};
use crate::error::{NmrIoError, Result};
use crate::virtual_fs::LoadedDataset;
use nmr_core::{AcquisitionParams, AxisMeta, Domain, NmrFormat, Signal1D, Signal2D};
use num_complex::Complex64;
use std::collections::BTreeMap;

/// Parse a JEOL Delta (`.jdf`) binary file in pure Rust.
///
/// Supports big-endian (`JEOL.NMR`) and little-endian (`RMN.LOEJ`) files, 1D & 2D
/// submatrix layouts, `float64`/`float32` data, and parameter block extraction.
pub fn load_jdf_bytes(filename: &str, bytes: &[u8]) -> Result<LoadedDataset> {
    if bytes.len() < 1400 {
        return Err(NmrIoError::Malformed {
            format: "JEOL .jdf",
            detail: "File too short for JEOL header".to_string(),
        });
    }

    let mut r = BinReader::new(bytes, Endian::Big, "JEOL .jdf");
    let ident = r.take(8)?;
    if ident != b"JEOL.NMR" && ident != b"RMN.LOEJ" {
        return Err(NmrIoError::Malformed {
            format: "JEOL .jdf",
            detail: "Missing JEOL.NMR/RMN.LOEJ magic bytes".to_string(),
        });
    }

    // Header scalars before multi-byte endian-dependent fields:
    let endian_flag = r.read_i8()?;
    let header_endian = if endian_flag == 1 {
        Endian::Little
    } else {
        Endian::Big
    };
    let _major_ver = r.read_u8()?;
    // Note: in JEOL spec, header numeric fields use the endianness in header_endian
    // except minor_version is read before/after — let's check both cleanly:
    let _minor_ver = r.read_u16()?;
    let _ndim_raw = r.read_u8()?;
    let dim_exist_byte = r.read_u8()?;
    let type_fmt_byte = r.read_u8()?;
    let data_type_code = type_fmt_byte >> 6; // 0 = f64, 1 = f32
    let data_format_code = type_fmt_byte & 0x3f; // 1 = 1D, 2 = 2D, 12 = small_2D
    let _instrument = r.read_i8()?;
    let _translate = r.take(8)?;
    let axis_types = r.take(8)?.to_vec();
    let _units = r.take(16)?; // 8 * (1 + 1) bytes
    let title = r.read_ascii_trimmed(124)?;
    let _axis_ranged = r.take(4)?;

    // Switch reader to the file's declared endianness for numeric arrays
    r.set_endian(header_endian);

    let mut data_points = [0u32; 8];
    for p in &mut data_points {
        *p = r.read_u32()?;
    }
    let mut offset_start = [0u32; 8];
    for p in &mut offset_start {
        *p = r.read_u32()?;
    }
    let mut offset_stop = [0u32; 8];
    for p in &mut offset_stop {
        *p = r.read_u32()?;
    }
    let mut _axis_start = [0f64; 8];
    for p in &mut _axis_start {
        *p = r.read_f64()?;
    }
    let mut _axis_stop = [0f64; 8];
    for p in &mut _axis_stop {
        *p = r.read_f64()?;
    }
    // creation_time (4B), revision_time (4B), node_name (16B), site (128B), author (128B), comment (128B)
    r.skip(4 + 4 + 16 + 128 + 128 + 128)?;
    // data_axis_titles (8 * 32B)
    let mut axis_titles = Vec::with_capacity(8);
    for _ in 0..8 {
        axis_titles.push(r.read_ascii_trimmed(32)?);
    }
    let mut base_freq = [0f64; 8];
    for p in &mut base_freq {
        *p = r.read_f64()?;
    }
    let mut _zero_point = [0f64; 8];
    for p in &mut _zero_point {
        *p = r.read_f64()?;
    }
    r.skip(8 + 3 + 1)?; // reversed (8B) + skip(3) + annotation_ok(1)
    let _history_used = r.read_u32()?;
    let _history_len = r.read_u32()?;
    let param_start = r.read_u32()? as usize;
    let _param_len = r.read_u32()?;
    r.skip(8 * 4 * 2)?; // list_start[8], list_length[8]
    let data_start = r.read_u32()? as usize;

    // Parse parameter table if valid
    let mut vendor_params = BTreeMap::new();
    if param_start > 0 && param_start + 16 <= bytes.len() {
        if let Ok(params) = parse_jeol_params(bytes, param_start, header_endian) {
            vendor_params = params;
        }
    }

    // Active dimensions
    let ndim = axis_types.iter().filter(|&&t| t != 0).count().max(1);
    let dim_names: Vec<char> = "xyzabcde"
        .chars()
        .enumerate()
        .filter(|(i, _)| (dim_exist_byte >> (7 - i)) & 1 == 1)
        .map(|(_, c)| c)
        .collect();
    let d0 = dim_names.first().copied().unwrap_or('x');

    let complex_dims = axis_types
        .iter()
        .take(ndim)
        .filter(|&&t| t == 3 || t == 4)
        .count()
        .min(1);
    let nsections = 1usize << complex_dims;

    let npts_total: usize = data_points.iter().take(ndim).map(|&p| (p as usize).max(1)).product();
    let elem_count = npts_total * nsections;

    r.set_pos(data_start)?;
    let mut raw_floats = Vec::with_capacity(elem_count);
    for _ in 0..elem_count {
        let v = if data_type_code == 1 {
            r.read_f32()? as f64
        } else {
            r.read_f64()?
        };
        raw_floats.push(v);
    }

    let sw_hz = vendor_params
        .get(&format!("{d0}_sweep"))
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(6000.0);
    let sf_mhz = vendor_params
        .get(&format!("{d0}_freq"))
        .and_then(|s| s.parse::<f64>().ok())
        .map(|hz| hz / 1e6)
        .filter(|&v| v > 0.0)
        .unwrap_or_else(|| if base_freq[0] > 0.0 { base_freq[0] } else { 500.0 });
    let car_ppm = vendor_params
        .get(&format!("{d0}_offset"))
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(5.0);
    let domain = match vendor_params.get(&format!("{d0}_domain")).map(|s| s.as_str()) {
        Some("Time") | None => Domain::Time,
        _ => Domain::Frequency,
    };

    let meta = AcquisitionParams {
        title: if title.is_empty() {
            filename.to_string()
        } else {
            title
        },
        format: NmrFormat::Jeol,
        axis: AxisMeta {
            nucleus: axis_titles
                .first()
                .filter(|s| !s.is_empty())
                .cloned()
                .unwrap_or_else(|| "1H".to_string()),
            spectrometer_freq_mhz: sf_mhz,
            spectral_width_hz: sw_hz,
            offset_ppm: None,
            carrier_ppm: car_ppm,
        },
        indirect_axis: None,
        group_delay: None,
        num_scans: vendor_params
            .get("scans")
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(1),
        vendor_params,
    };

    if ndim >= 2 && (data_format_code == 2 || data_format_code == 12) {
        let cols = (data_points[0] as usize).max(1);
        let rows = (data_points[1] as usize).max(1);
        let edge = if data_format_code == 12 { 4 } else { 32 };
        let sec_len = rows * cols;
        let re_reordered = reorder_submatrix_2d(&raw_floats[..sec_len.min(raw_floats.len())], rows, cols, edge);
        let im_reordered = if nsections >= 2 && raw_floats.len() >= 2 * sec_len {
            reorder_submatrix_2d(&raw_floats[sec_len..2 * sec_len], rows, cols, edge)
        } else {
            vec![0.0; rows * cols]
        };
        let data: Vec<Complex64> = re_reordered
            .into_iter()
            .zip(im_reordered)
            .map(|(re, im)| Complex64::new(re, im))
            .collect();
        Ok(LoadedDataset::D2(Signal2D::new(
            rows,
            cols,
            data,
            domain,
            Domain::Time,
            meta,
        )))
    } else {
        let len = (data_points[0] as usize).min(raw_floats.len());
        let valid_start = offset_start[0] as usize;
        let valid_stop = (offset_stop[0] as usize).min(len.saturating_sub(1));
        let slice_end = if valid_stop > valid_start {
            valid_stop + 1
        } else {
            len
        };
        let mut data = Vec::with_capacity(slice_end.saturating_sub(valid_start));
        for idx in valid_start..slice_end {
            let re = raw_floats[idx];
            let im = if nsections >= 2 && len + idx < raw_floats.len() {
                raw_floats[len + idx]
            } else {
                0.0
            };
            data.push(Complex64::new(re, im));
        }
        Ok(LoadedDataset::D1(Signal1D::new(data, domain, meta)))
    }
}

fn parse_jeol_params(
    bytes: &[u8],
    param_start: usize,
    endian: Endian,
) -> Result<BTreeMap<String, String>> {
    let mut r = BinReader::new(bytes, endian, "JEOL params");
    r.set_pos(param_start)?;
    let _param_size = r.read_u32()?;
    let _low_idx = r.read_u32()?;
    let high_idx = r.read_u32()? as usize;
    let _total_size = r.read_u32()?;

    let mut map = BTreeMap::new();
    for _ in 0..high_idx.min(2048) {
        if r.pos() + 64 > bytes.len() {
            break;
        }
        r.skip(4)?; // _class
        let unit_scaler = r.read_i16()?;
        r.skip(10)?; // unit[5] * 2B
        r.skip(16)?; // value bytes preview
        let val_type = r.read_i32()?;
        // rewind 20 bytes to read value properly
        r.set_pos(r.pos() - 20)?;
        let scale = 10.0_f64.powi(unit_scaler as i32);
        let val_str = match val_type {
            0 => r.read_ascii_trimmed(16)?,
            1 | 4 => {
                let v = r.read_i32()?;
                r.skip(12)?;
                ((v as f64) * scale).to_string()
            }
            2 => {
                let v = r.read_f64()?;
                r.skip(8)?;
                (v * scale).to_string()
            }
            3 => {
                let re = r.read_f64()?;
                let im = r.read_f64()?;
                format!("({},{})", re * scale, im * scale)
            }
            _ => {
                r.skip(16)?;
                String::new()
            }
        };
        r.skip(4)?; // val_type again
        let name = r.read_ascii_trimmed(28)?;
        let key = name.to_ascii_lowercase().replace(' ', "");
        if !key.is_empty() {
            map.insert(key, val_str);
        }
    }
    Ok(map)
}

/// Reorder 2D tiled JEOL submatrix data (`edge × edge` blocks) into row-major `rows × cols`.
fn reorder_submatrix_2d(src: &[f64], rows: usize, cols: usize, edge: usize) -> Vec<f64> {
    let mut out = vec![0.0; rows * cols];
    let blocks_r = rows.div_ceil(edge);
    let blocks_c = cols.div_ceil(edge);
    let mut idx = 0usize;
    for br in 0..blocks_r {
        for bc in 0..blocks_c {
            for r in 0..edge {
                for c in 0..edge {
                    if idx >= src.len() {
                        return out;
                    }
                    let gr = br * edge + r;
                    let gc = bc * edge + c;
                    if gr < rows && gc < cols {
                        out[gr * cols + gc] = src[idx];
                    }
                    idx += 1;
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loads_minimal_synthetic_jdf_1d() {
        let mut buf = vec![0u8; 1600];
        buf[0..8].copy_from_slice(b"JEOL.NMR");
        buf[8] = 1; // little endian
        buf[9] = 1; // major
        buf[12] = 1; // ndim
        buf[13] = 0b1000_0000; // x exists
        buf[14] = (0 << 6) | 1; // float64, one_d
        buf[24] = 3; // x axis complex
        // title at offset 48 (after 8+1+1+2+1+1+1+1+8+8+16 = 48)
        buf[48..57].copy_from_slice(b"Ethanol1H");
        // data_points[0] at offset 48 + 124 + 4 = 176
        buf[176..180].copy_from_slice(&4u32.to_le_bytes());
        // offset_stop[0] at 176 + 32 + 32 = 240
        buf[240..244].copy_from_slice(&3u32.to_le_bytes());
        // base_freq[0] at 176 + 32*3 + 64*2 + 4+4+16+128+128+128 + 256 = 176 + 96 + 128 + 408 + 256 = 1064
        buf[1064..1072].copy_from_slice(&500.0_f64.to_le_bytes());
        // data_start offset: 1064 + 64(base_freq) + 64(zero_point) + 12(rev/ann) + 8(hist) + 8(param) + 64(list) = 1284
        let data_start: usize = 1400;
        buf[1284..1288].copy_from_slice(&(data_start as u32).to_le_bytes());
        // Write 4 real + 4 imag f64 values
        for i in 0..4 {
            let off_re = data_start + i * 8;
            let off_im = data_start + (4 + i) * 8;
            buf[off_re..off_re + 8].copy_from_slice(&((i as f64) + 1.0).to_le_bytes());
            buf[off_im..off_im + 8].copy_from_slice(&(-((i as f64) + 0.5)).to_le_bytes());
        }

        let ds = load_jdf_bytes("sample.jdf", &buf).unwrap().into_1d().unwrap();
        assert_eq!(ds.len(), 4);
        assert_eq!(ds.meta.title, "Ethanol1H");
        assert!((ds.data[1].re - 2.0).abs() < 1e-9);
        assert!((ds.data[1].im + 1.5).abs() < 1e-9);
    }
}
