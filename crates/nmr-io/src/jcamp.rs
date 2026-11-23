use crate::error::{NmrIoError, Result};
use nmr_core::{AcquisitionParams, AxisMeta, Domain, NmrFormat, Signal1D};
use num_complex::Complex64;
use std::collections::BTreeMap;

/// Parse a JCAMP-DX (`.dx` / `.jdx`) text buffer, supporting AFFN numbers,
/// ASDF compressed encoding (SQZ, DIF, DUP), and multi-page `NTUPLES` (e.g. Real + Imaginary).
pub fn load_jcamp_str(filename: &str, text: &str) -> Result<Signal1D> {
    let mut headers = BTreeMap::<String, String>::new();
    let mut pages: Vec<(String, Vec<f64>)> = Vec::new();
    let mut current_page_symbol = String::from("Y");
    let mut current_y: Vec<f64> = Vec::new();
    let mut in_xydata = false;
    let mut page_xfactor: Option<f64> = None;
    let mut page_yfactor: Option<f64> = None;
    let mut page_firstx: Option<f64> = None;
    let mut page_lastx: Option<f64> = None;

    for raw_line in text.lines() {
        let line = raw_line.split("$$").next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        if let Some(rest) = line.strip_prefix("##") {
            if in_xydata && !current_y.is_empty() {
                let yf = page_yfactor
                    .or_else(|| headers.get("YFACTOR").and_then(|s| s.parse().ok()))
                    .unwrap_or(1.0);
                let scaled: Vec<f64> = current_y.iter().map(|v| v * yf).collect();
                pages.push((current_page_symbol.clone(), scaled));
                current_y.clear();
            }
            in_xydata = false;

            if rest.eq_ignore_ascii_case("END=") {
                break;
            }

            if let Some((k_raw, v_raw)) = rest.split_once('=') {
                let key = k_raw.trim().to_ascii_uppercase();
                let val = v_raw.trim();
                if key == "XYDATA" || key == "DATA TABLE" {
                    in_xydata = true;
                    continue;
                }
                if key == "PAGE" {
                    // e.g. ##PAGE= N=1 or ##PAGE= SYMBOL= R
                    if let Some(idx) = val.to_ascii_uppercase().find("SYMBOL=") {
                        current_page_symbol = val[idx + 7..].split(',').next().unwrap_or("Y").trim().to_string();
                    } else if let Some(idx) = val.to_ascii_uppercase().find("N=") {
                        current_page_symbol = format!("PAGE_{}", val[idx + 2..].trim());
                    }
                    continue;
                }
                match key.as_str() {
                    "XFACTOR" => page_xfactor = val.parse().ok(),
                    "YFACTOR" => page_yfactor = val.parse().ok(),
                    "FIRSTX" => page_firstx = val.parse().ok(),
                    "LASTX" => page_lastx = val.parse().ok(),
                    _ => {}
                }
                headers.entry(key).or_insert_with(|| val.to_string());
            }
        } else if in_xydata {
            let mut decoded = parse_asdf_xy_line(line)?;
            current_y.append(&mut decoded);
        }
    }

    if !current_y.is_empty() {
        let yf = page_yfactor
            .or_else(|| headers.get("YFACTOR").and_then(|s| s.parse().ok()))
            .unwrap_or(1.0);
        let scaled: Vec<f64> = current_y.iter().map(|v| v * yf).collect();
        pages.push((current_page_symbol, scaled));
    }

    if pages.is_empty() {
        return Err(NmrIoError::Malformed {
            format: "JCAMP-DX",
            detail: "No XYDATA or DATA TABLE samples found".to_string(),
        });
    }

    let data_type = headers
        .get("DATA TYPE")
        .map(|s| s.to_ascii_uppercase())
        .unwrap_or_default();
    let x_units = headers
        .get("XUNITS")
        .map(|s| s.to_ascii_uppercase())
        .unwrap_or_default();
    let domain = if data_type.contains("FID") || x_units.contains("SEC") {
        Domain::Time
    } else {
        Domain::Frequency
    };

    let n = pages[0].1.len();
    let data: Vec<Complex64> = if pages.len() >= 2 {
        // Combine Real + Imaginary pages (common in NTUPLES FID/Spectrum)
        (0..n)
            .map(|i| {
                let re = pages[0].1[i];
                let im = pages[1].1.get(i).copied().unwrap_or(0.0);
                Complex64::new(re, im)
            })
            .collect()
    } else {
        pages[0].1.iter().map(|&r| Complex64::new(r, 0.0)).collect()
    };

    let sf = headers
        .get(".OBSERVE FREQUENCY")
        .or_else(|| headers.get("$BF1"))
        .or_else(|| headers.get("$SF"))
        .and_then(|s| s.split_whitespace().next())
        .and_then(|s| s.parse::<f64>().ok())
        .filter(|&v| v > 0.0)
        .unwrap_or(500.0);

    let first_x = page_firstx
        .or_else(|| headers.get("FIRSTX").and_then(|s| s.parse().ok()))
        .unwrap_or(10.0);
    let last_x = page_lastx
        .or_else(|| headers.get("LASTX").and_then(|s| s.parse().ok()))
        .unwrap_or(0.0);
    let _xf = page_xfactor.unwrap_or(1.0);

    let (sw_hz, offset_ppm, carrier_ppm) = if domain == Domain::Frequency {
        if x_units.contains("HZ") {
            let sw = (first_x - last_x).abs().max(1.0);
            let left = first_x.max(last_x) / sf;
            (sw, Some(left), left - 0.5 * (sw / sf))
        } else {
            // PPM default
            let sw_ppm = (first_x - last_x).abs().max(1e-6);
            let sw = sw_ppm * sf;
            let left = first_x.max(last_x);
            (sw, Some(left), left - 0.5 * sw_ppm)
        }
    } else {
        let dt = ((last_x - first_x).abs() / (n.saturating_sub(1).max(1)) as f64).max(1e-9);
        let sw = 1.0 / dt;
        (sw, None, 5.0)
    };

    let title = headers
        .get("TITLE")
        .cloned()
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| filename.to_string());
    let nucleus = headers
        .get(".OBSERVE NUCLEUS")
        .map(|s| s.trim_start_matches('^').to_string())
        .unwrap_or_else(|| "1H".to_string());

    let meta = AcquisitionParams {
        title,
        format: NmrFormat::JcampDx,
        axis: AxisMeta {
            nucleus,
            spectrometer_freq_mhz: sf,
            spectral_width_hz: sw_hz,
            offset_ppm,
            carrier_ppm,
        },
        indirect_axis: None,
        group_delay: None,
        num_scans: 1,
        vendor_params: headers,
    };

    Ok(Signal1D::new(data, domain, meta))
}

/// Decode one JCAMP-DX `(X++(Y..Y))` line supporting AFFN numbers and ASDF SQZ/DIF/DUP compression.
///
/// Reference encoding rules per IUPAC JCAMP-DX specification:
/// - Leading number on the line is the X coordinate (ignored because Y values are uniformly spaced).
/// - SQZ: `@`=`+0`, `A..I`=`+1..+9`, `a..i`=`-1..-9`
/// - DIF: `%`=`+0`, `J..R`=`+1..+9`, `j..r`=`-1..-9`
/// - DUP: `S..Z`=`1..8`, `s`=`9` (multi-digit counts e.g. `s0` = 90)
fn parse_asdf_xy_line(line: &str) -> Result<Vec<f64>> {
    // Tokenize skipping the leading X coordinate.
    let bytes = line.as_bytes();
    let mut pos = 0usize;

    // Skip leading whitespace and the initial X value (digits, sign, '.', 'E'/'e').
    while pos < bytes.len() && bytes[pos].is_ascii_whitespace() {
        pos += 1;
    }
    while pos < bytes.len()
        && (bytes[pos].is_ascii_digit()
            || matches!(bytes[pos], b'+' | b'-' | b'.' | b'e' | b'E'))
    {
        pos += 1;
    }

    let mut y_vals: Vec<f64> = Vec::new();
    let mut last_dif: Option<i64> = None;

    // Walk remaining characters; handles both space/comma-separated AFFN and packed ASDF.
    let rest = &line[pos..];
    if !rest.bytes().any(is_asdf_pseudo_digit) {
        // Pure AFFN numbers separated by spaces/commas or +/-
        for tok in rest.split(|c: char| c.is_ascii_whitespace() || c == ',' || c == ';') {
            let t = tok.trim();
            if !t.is_empty() {
                if let Ok(v) = t.parse::<f64>() {
                    y_vals.push(v);
                }
            }
        }
        return Ok(y_vals);
    }

    // ASDF decoder
    let chars: Vec<char> = rest.chars().collect();
    let mut i = 0usize;
    let mut last_val: f64 = 0.0;

    while i < chars.len() {
        let ch = chars[i];
        if ch.is_ascii_whitespace() || ch == ',' {
            i += 1;
            continue;
        }
        if let Some(dup_first) = dup_digit(ch) {
            // DUP repeat count
            let mut count_str = dup_first.to_string();
            i += 1;
            while i < chars.len() && chars[i].is_ascii_digit() {
                count_str.push(chars[i]);
                i += 1;
            }
            let count: usize = count_str.parse().unwrap_or(1);
            let repeat = count.saturating_sub(1);
            for _ in 0..repeat {
                if let Some(d) = last_dif {
                    last_val += d as f64;
                }
                y_vals.push(last_val);
            }
            continue;
        }

        if let Some((is_dif, first_digit)) = sqz_or_dif_digit(ch) {
            i += 1;
            let mut num_str = first_digit;
            while i < chars.len() && (chars[i].is_ascii_digit() || chars[i] == '.') {
                num_str.push(chars[i]);
                i += 1;
            }
            let parsed: f64 = num_str.parse().map_err(|_| NmrIoError::Malformed {
                format: "JCAMP-DX",
                detail: format!("Invalid ASDF token '{num_str}'"),
            })?;
            if is_dif {
                let d = parsed.round() as i64;
                last_val += d as f64;
                last_dif = Some(d);
            } else {
                last_val = parsed;
                last_dif = None;
            }
            y_vals.push(last_val);
            continue;
        }

        // Plain signed/unsigned number inside ASDF line
        if ch.is_ascii_digit() || matches!(ch, '+' | '-' | '.') {
            let start = i;
            i += 1;
            while i < chars.len()
                && (chars[i].is_ascii_digit() || matches!(chars[i], '.' | 'e' | 'E'))
            {
                if matches!(chars[i], 'e' | 'E') && i + 1 < chars.len() && matches!(chars[i + 1], '+' | '-') {
                    i += 2;
                } else {
                    i += 1;
                }
            }
            let tok: String = chars[start..i].iter().collect();
            let parsed: f64 = tok.parse().map_err(|_| NmrIoError::Malformed {
                format: "JCAMP-DX",
                detail: format!("Invalid AFFN number '{tok}'"),
            })?;
            last_val = parsed;
            last_dif = None;
            y_vals.push(last_val);
            continue;
        }

        i += 1;
    }

    Ok(y_vals)
}

fn is_asdf_pseudo_digit(b: u8) -> bool {
    matches!(b, b'@' | b'%' | b'A'..=b'Z' | b'a'..=b'i' | b'j'..=b'r' | b's')
}

fn sqz_or_dif_digit(c: char) -> Option<(bool, String)> {
    match c {
        '@' => Some((false, "0".into())),
        'A'..='I' => Some((false, ((c as u8 - b'A') + 1).to_string())),
        'a'..='i' => Some((false, format!("-{}", (c as u8 - b'a') + 1))),
        '%' => Some((true, "0".into())),
        'J'..='R' => Some((true, ((c as u8 - b'J') + 1).to_string())),
        'j'..='r' => Some((true, format!("-{}", (c as u8 - b'j') + 1))),
        _ => None,
    }
}

fn dup_digit(c: char) -> Option<u8> {
    match c {
        'S'..='Z' => Some((c as u8 - b'S') + 1),
        's' => Some(9),
        _ => None,
    }
}

/// Export a [`Signal1D`] to a standard AFFN JCAMP-DX string (supports Real + Imaginary `NTUPLES`).
pub fn write_jcamp_dx(sig: &Signal1D) -> String {
    use std::fmt::Write;
    let mut out = String::new();
    let n = sig.len();
    let (data_type, x_units, first_x, last_x) = match sig.domain {
        Domain::Frequency => (
            "NMR SPECTRUM",
            "PPM",
            sig.meta.axis.left_ppm(),
            sig.meta.axis.right_ppm(),
        ),
        Domain::Time => (
            "NMR FID",
            "SECONDS",
            0.0,
            sig.meta.axis.index_to_time_s(n.saturating_sub(1)),
        ),
    };
    let delta_x = if n > 1 {
        (last_x - first_x) / (n - 1) as f64
    } else {
        0.0
    };

    let _ = writeln!(out, "##TITLE= {}", sig.meta.title);
    let _ = writeln!(out, "##JCAMP-DX= 5.00");
    let _ = writeln!(out, "##DATA TYPE= {data_type}");
    let _ = writeln!(out, "##DATA CLASS= NTUPLES");
    let _ = writeln!(out, "##.OBSERVE NUCLEUS= ^{}", sig.meta.axis.nucleus);
    let _ = writeln!(out, "##.OBSERVE FREQUENCY= {:.6}", sig.meta.axis.spectrometer_freq_mhz);
    let _ = writeln!(out, "##NTUPLES= {data_type}");
    let _ = writeln!(out, "##VAR_NAME=  X, REAL, IMAG");
    let _ = writeln!(out, "##SYMBOL=    X, R, I");
    let _ = writeln!(out, "##VAR_DIM=   {n}, {n}, {n}");
    let _ = writeln!(out, "##UNITS=     {x_units}, ARBITRARY UNITS, ARBITRARY UNITS");
    let _ = writeln!(out, "##FIRST=     {first_x:.6}, 0, 0");
    let _ = writeln!(out, "##LAST=      {last_x:.6}, 0, 0");
    let _ = writeln!(out, "##XFACTOR=   1.0");
    let _ = writeln!(out, "##YFACTOR=   1.0");
    let _ = writeln!(out, "##FIRSTX=    {first_x:.6}");
    let _ = writeln!(out, "##LASTX=     {last_x:.6}");
    let _ = writeln!(out, "##NPOINTS=   {n}");

    for (page_idx, sym, getter) in [
        (1, "R", Box::new(|c: Complex64| c.re) as Box<dyn Fn(Complex64) -> f64>),
        (2, "I", Box::new(|c: Complex64| c.im) as Box<dyn Fn(Complex64) -> f64>),
    ] {
        let _ = writeln!(out, "##PAGE= N={page_idx} SYMBOL={sym}");
        let _ = writeln!(out, "##DATA TABLE= (X++(Y..Y)), XYDATA");
        for chunk_start in (0..n).step_by(8) {
            let x = first_x + chunk_start as f64 * delta_x;
            let _ = write!(out, "{x:.6}");
            for &c in &sig.data[chunk_start..(chunk_start + 8).min(n)] {
                let _ = write!(out, " {:.8}", getter(c));
            }
            out.push('\n');
        }
    }
    let _ = writeln!(out, "##END NTUPLES= {data_type}");
    let _ = writeln!(out, "##END=");
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_asdf_sqz_dif_dup_and_ntuples() {
        // First line: X=10.0, Y=100 (@=0 -> A00 = +100), then DIF J0 (+10 -> 110), DUP T (repeat 2 -> 120), SQZ b0 (-20)
        let dx = "\
##TITLE= Test ASDF
##JCAMP-DX= 5.00
##DATA TYPE= NMR SPECTRUM
##XUNITS= PPM
##.OBSERVE FREQUENCY= 400.0
##FIRSTX= 10.0
##LASTX= 2.0
##YFACTOR= 0.5
##PAGE= N=1 SYMBOL=R
##XYDATA= (X++(Y..Y))
10.0 A00J0Tb0
##PAGE= N=2 SYMBOL=I
##XYDATA= (X++(Y..Y))
10.0 10.0 20.0 30.0 40.0
##END=
";
        let sig = load_jcamp_str("test.jdx", dx).unwrap();
        assert_eq!(sig.len(), 4);
        assert_eq!(sig.domain, Domain::Frequency);
        assert!((sig.data[0].re - 50.0).abs() < 1e-9); // 100 * 0.5
        assert!((sig.data[1].re - 55.0).abs() < 1e-9); // 110 * 0.5
        assert!((sig.data[2].re - 60.0).abs() < 1e-9); // 120 * 0.5
        assert!((sig.data[3].re + 10.0).abs() < 1e-9); // -20 * 0.5
        assert!((sig.data[1].im - 10.0).abs() < 1e-9); // 20 * 0.5
    }

    #[test]
    fn roundtrip_jcamp_dx() {
        let mut meta = AcquisitionParams::default();
        meta.title = "Roundtrip".into();
        let sig = Signal1D::new(
            vec![Complex64::new(1.5, -0.5), Complex64::new(2.25, 0.75)],
            Domain::Frequency,
            meta,
        );
        let text = write_jcamp_dx(&sig);
        let loaded = load_jcamp_str("out.jdx", &text).unwrap();
        assert_eq!(loaded.len(), 2);
        assert!((loaded.data[0].re - 1.5).abs() < 1e-6);
        assert!((loaded.data[1].im - 0.75).abs() < 1e-6);
    }
}
