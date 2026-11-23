use crate::error::{NmrIoError, Result};

/// Endianness helper for multi-vendor NMR binary formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Endian {
    Little,
    Big,
}

/// Cursor-based binary reader supporting dynamic endianness.
pub struct BinReader<'a> {
    data: &'a [u8],
    pos: usize,
    endian: Endian,
    format: &'static str,
}

impl<'a> BinReader<'a> {
    pub fn new(data: &'a [u8], endian: Endian, format: &'static str) -> Self {
        Self {
            data,
            pos: 0,
            endian,
            format,
        }
    }

    pub fn set_endian(&mut self, endian: Endian) {
        self.endian = endian;
    }

    pub fn pos(&self) -> usize {
        self.pos
    }

    pub fn set_pos(&mut self, pos: usize) -> Result<()> {
        if pos > self.data.len() {
            return Err(NmrIoError::Malformed {
                format: self.format,
                detail: format!("Seek offset {pos} exceeds buffer length {}", self.data.len()),
            });
        }
        self.pos = pos;
        Ok(())
    }

    pub fn skip(&mut self, n: usize) -> Result<()> {
        self.set_pos(self.pos.saturating_add(n))
    }

    pub fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        let end = self.pos.checked_add(n).ok_or_else(|| NmrIoError::Malformed {
            format: self.format,
            detail: "Integer overflow reading bytes".to_string(),
        })?;
        if end > self.data.len() {
            return Err(NmrIoError::Malformed {
                format: self.format,
                detail: format!(
                    "Unexpected EOF at offset {} (needed {n} bytes, have {})",
                    self.pos,
                    self.data.len() - self.pos
                ),
            });
        }
        let slice = &self.data[self.pos..end];
        self.pos = end;
        Ok(slice)
    }

    pub fn read_u8(&mut self) -> Result<u8> {
        Ok(self.take(1)?[0])
    }

    pub fn read_i8(&mut self) -> Result<i8> {
        Ok(self.read_u8()? as i8)
    }

    pub fn read_i16(&mut self) -> Result<i16> {
        let b: [u8; 2] = self.take(2)?.try_into().unwrap();
        Ok(match self.endian {
            Endian::Little => i16::from_le_bytes(b),
            Endian::Big => i16::from_be_bytes(b),
        })
    }

    pub fn read_u16(&mut self) -> Result<u16> {
        let b: [u8; 2] = self.take(2)?.try_into().unwrap();
        Ok(match self.endian {
            Endian::Little => u16::from_le_bytes(b),
            Endian::Big => u16::from_be_bytes(b),
        })
    }

    pub fn read_i32(&mut self) -> Result<i32> {
        let b: [u8; 4] = self.take(4)?.try_into().unwrap();
        Ok(match self.endian {
            Endian::Little => i32::from_le_bytes(b),
            Endian::Big => i32::from_be_bytes(b),
        })
    }

    pub fn read_u32(&mut self) -> Result<u32> {
        let b: [u8; 4] = self.take(4)?.try_into().unwrap();
        Ok(match self.endian {
            Endian::Little => u32::from_le_bytes(b),
            Endian::Big => u32::from_be_bytes(b),
        })
    }

    pub fn read_f32(&mut self) -> Result<f32> {
        let b: [u8; 4] = self.take(4)?.try_into().unwrap();
        Ok(match self.endian {
            Endian::Little => f32::from_le_bytes(b),
            Endian::Big => f32::from_be_bytes(b),
        })
    }

    pub fn read_f64(&mut self) -> Result<f64> {
        let b: [u8; 8] = self.take(8)?.try_into().unwrap();
        Ok(match self.endian {
            Endian::Little => f64::from_le_bytes(b),
            Endian::Big => f64::from_be_bytes(b),
        })
    }

    pub fn read_ascii_trimmed(&mut self, n: usize) -> Result<String> {
        let bytes = self.take(n)?;
        let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
        Ok(String::from_utf8_lossy(&bytes[..end]).trim().to_string())
    }
}
