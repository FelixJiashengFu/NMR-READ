use crate::error::{NmrIoError, Result};
use base64::{Engine as _, engine::general_purpose::STANDARD as BASE64};
use nmr_core::{AcquisitionParams, AxisMeta, Domain, NmrFormat, Signal1D};
use num_complex::Complex64;
use quick_xml::Reader;
use quick_xml::events::Event;
use std::collections::BTreeMap;

/// Parse an open standard `nmrML` XML document containing `<fidData>` or `<spectrumDataArray>`.
pub fn load_nmrml_str(filename: &str, xml: &str) -> Result<Signal1D> {
    let mut reader = Reader::from_str(xml);
    let mut buf = Vec::new();
    let mut in_fid_data = false;
    let mut in_spec_data = false;
    let mut byte_format = String::from("complex128");
    let mut encoded_text = String::new();
    let mut is_frequency = false;

    let mut sf_mhz = 500.0_f64;
    let mut sw_hz = 6000.0_f64;
    let mut nucleus = String::from("1H");
    let mut vendor_params = BTreeMap::new();

    loop {
        match reader.read_event_into(&mut buf) {
            Ok(Event::Start(e)) | Ok(Event::Empty(e)) => {
                let qname = e.name();
                let tag = qname.as_ref();
                let local = tag.rsplit(':').next().unwrap_or(tag);
                if local == "fidData" {
                    in_fid_data = true;
                    is_frequency = false;
                    encoded_text.clear();
                    for attr in e.attributes().flatten() {
                        if attr.key.as_ref() == "byteFormat" {
                            byte_format = attr.value.to_string();
                        }
                    }
                } else if local == "spectrumDataArray" {
                    in_spec_data = true;
                    is_frequency = true;
                    encoded_text.clear();
                    for attr in e.attributes().flatten() {
                        if attr.key.as_ref() == "byteFormat" {
                            byte_format = attr.value.to_string();
                        }
                    }
                } else if local == "effectiveExcitationFrequency" || local == "irradiationFrequency" {
                    for attr in e.attributes().flatten() {
                        if attr.key.as_ref() == "value" {
                            if let Ok(v) = attr.value.parse::<f64>() {
                                // Convert Hz -> MHz if large
                                sf_mhz = if v > 1e5 { v / 1e6 } else { v };
                                vendor_params.insert(local.to_string(), v.to_string());
                            }
                        }
                    }
                } else if local == "sweepWidth" {
                    for attr in e.attributes().flatten() {
                        if attr.key.as_ref() == "value" {
                            if let Ok(v) = attr.value.parse::<f64>() {
                                sw_hz = v;
                                vendor_params.insert("sweepWidth".to_string(), v.to_string());
                            }
                        }
                    }
                } else if local == "acquisitionNucleus" {
                    for attr in e.attributes().flatten() {
                        if attr.key.as_ref() == "name" {
                            nucleus = attr.value.to_string();
                        }
                    }
                }
            }
            Ok(Event::Text(t)) => {
                if in_fid_data || in_spec_data {
                    let s = t.into_inner();
                    encoded_text.push_str(s.trim());
                }
            }
            Ok(Event::End(e)) => {
                let qname = e.name();
                let tag = qname.as_ref();
                let local = tag.rsplit(':').next().unwrap_or(tag);
                if local == "fidData" || local == "spectrumDataArray" {
                    in_fid_data = false;
                    in_spec_data = false;
                    if !encoded_text.is_empty() {
                        break;
                    }
                }
            }
            Ok(Event::Eof) => break,
            Err(e) => return Err(NmrIoError::Xml(e.to_string())),
            _ => {}
        }
        buf.clear();
    }

    if encoded_text.is_empty() {
        return Err(NmrIoError::Malformed {
            format: "nmrML",
            detail: "No <fidData> or <spectrumDataArray> payload found".to_string(),
        });
    }

    let clean_b64: String = encoded_text.chars().filter(|c| !c.is_whitespace()).collect();
    let raw = BASE64.decode(clean_b64).map_err(|e| NmrIoError::Malformed {
        format: "nmrML",
        detail: format!("Invalid base64 payload: {e}"),
    })?;

    let data: Vec<Complex64> = if byte_format.contains("64") && !byte_format.contains("complex128") {
        // Real float64 array
        raw.chunks_exact(8)
            .map(|c| Complex64::new(f64::from_le_bytes(c.try_into().unwrap()), 0.0))
            .collect()
    } else if byte_format.contains("32") || byte_format.contains("complex64") {
        // Interleaved f32 complex pairs
        raw.chunks_exact(8)
            .map(|c| {
                let re = f32::from_le_bytes(c[0..4].try_into().unwrap()) as f64;
                let im = f32::from_le_bytes(c[4..8].try_into().unwrap()) as f64;
                Complex64::new(re, im)
            })
            .collect()
    } else {
        // Default: complex128 (interleaved little-endian f64 Re, Im)
        raw.chunks_exact(16)
            .map(|c| {
                let re = f64::from_le_bytes(c[0..8].try_into().unwrap());
                let im = f64::from_le_bytes(c[8..16].try_into().unwrap());
                Complex64::new(re, im)
            })
            .collect()
    };

    let meta = AcquisitionParams {
        title: filename.to_string(),
        format: NmrFormat::NmrMl,
        axis: AxisMeta {
            nucleus,
            spectrometer_freq_mhz: sf_mhz,
            spectral_width_hz: sw_hz,
            offset_ppm: None,
            carrier_ppm: 5.0,
        },
        indirect_axis: None,
        group_delay: None,
        num_scans: 1,
        vendor_params,
    };

    Ok(Signal1D::new(
        data,
        if is_frequency {
            Domain::Frequency
        } else {
            Domain::Time
        },
        meta,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loads_nmrml_complex128_fid() {
        let mut raw = Vec::new();
        for (re, im) in [(1.0_f64, -2.0_f64), (3.0, 4.0)] {
            raw.extend_from_slice(&re.to_le_bytes());
            raw.extend_from_slice(&im.to_le_bytes());
        }
        let b64 = BASE64.encode(raw);
        let xml = format!(
            r#"<?xml version="1.0"?>
<nmrML>
  <acquisition>
    <acquisitionParameterSet>
      <DirectDimensionParameterSet>
        <acquisitionNucleus name="1H"/>
        <effectiveExcitationFrequency value="500.13"/>
        <sweepWidth value="6000.0"/>
      </DirectDimensionParameterSet>
    </acquisitionParameterSet>
    <fidData byteFormat="complex128">{b64}</fidData>
  </acquisition>
</nmrML>"#
        );
        let sig = load_nmrml_str("test.nmrML", &xml).unwrap();
        assert_eq!(sig.len(), 2);
        assert_eq!(sig.domain, Domain::Time);
        assert!((sig.meta.axis.spectrometer_freq_mhz - 500.13).abs() < 1e-6);
        assert!((sig.data[1].re - 3.0).abs() < 1e-9);
        assert!((sig.data[1].im - 4.0).abs() < 1e-9);
    }
}
