//! Enlace serie con buffering y encuadrado de tramas para comunicación PC ↔ ESP32.

use crate::protocol::{Frame, CRC_SIZE, HEADER_SIZE, SYNC_WORD};
use std::io::{Read, Write};
use std::time::{Duration, Instant};

pub struct SerialLink {
    port: Box<dyn serialport::SerialPort>,
    rx_buf: Vec<u8>,
}

impl SerialLink {
    /// Abre un puerto serie con la configuración especificada.
    pub fn open(path: &str, baudrate: u32) -> Result<Self, String> {
        let port = serialport::new(path, baudrate)
            .timeout(Duration::from_millis(10))
            .data_bits(serialport::DataBits::Eight)
            .stop_bits(serialport::StopBits::One)
            .parity(serialport::Parity::None)
            .flow_control(serialport::FlowControl::None)
            .open()
            .map_err(|e| format!("No se pudo abrir {}: {}", path, e))?;

        Ok(SerialLink {
            port,
            rx_buf: Vec::with_capacity(4096),
        })
    }

    /// Escribe una trama por el enlace serie.
    pub fn write_frame(&mut self, frame: &Frame) -> Result<(), String> {
        let bytes = frame.to_bytes();
        self.port
            .write_all(&bytes)
            .map_err(|e| format!("Error escritura: {}", e))?;
        self.port
            .flush()
            .map_err(|e| format!("Error flush: {}", e))?;
        Ok(())
    }

    /// Lee una trama completa del enlace serie con timeout.
    pub fn read_frame(&mut self, timeout_ms: u64) -> Result<Option<Frame>, String> {
        let deadline = Instant::now() + Duration::from_millis(timeout_ms);
        let mut tmp = [0u8; 256];

        loop {
            // Intentar decodificar lo que ya tenemos en buffer
            if let Some(frame) = self.try_decode()? {
                return Ok(Some(frame));
            }

            if Instant::now() > deadline {
                return Ok(None); // timeout sin trama completa
            }

            match self.port.read(&mut tmp) {
                Ok(n) if n > 0 => {
                    self.rx_buf.extend_from_slice(&tmp[..n]);
                    // Limitar tamaño del buffer para evitar crecimiento ilimitado
                    if self.rx_buf.len() > 64 * 1024 {
                        self.rx_buf.drain(..32 * 1024);
                    }
                }
                Ok(_) => {
                    // 0 bytes leídos, esperar un poco
                    std::thread::sleep(Duration::from_millis(1));
                }
                Err(ref e) if e.kind() == std::io::ErrorKind::TimedOut => {
                    std::thread::sleep(Duration::from_millis(1));
                }
                Err(e) => return Err(format!("Error lectura: {}", e)),
            }
        }
    }

    /// Intenta extraer una trama del buffer de recepción.
    fn try_decode(&mut self) -> Result<Option<Frame>, String> {
        // Buscar sync word
        let sync_pos = self
            .rx_buf
            .windows(2)
            .position(|w| w[0] == SYNC_WORD[0] && w[1] == SYNC_WORD[1]);

        let pos = match sync_pos {
            Some(p) => p,
            None => {
                // Si no encontramos sync word y el buffer es grande, descartamos
                if self.rx_buf.len() > 2 {
                    self.rx_buf.drain(..self.rx_buf.len() - 1);
                }
                return Ok(None);
            }
        };

        // Descartar bytes antes de la sync word
        if pos > 0 {
            self.rx_buf.drain(..pos);
        }

        // Necesitamos al menos header para saber la longitud
        if self.rx_buf.len() < HEADER_SIZE {
            return Ok(None);
        }

        let payload_len = u16::from_le_bytes([self.rx_buf[2], self.rx_buf[3]]) as usize;
        let total_len = HEADER_SIZE + payload_len + CRC_SIZE;

        if self.rx_buf.len() < total_len {
            return Ok(None);
        }

        // Extraer trama completa
        let frame_bytes: Vec<u8> = self.rx_buf.drain(..total_len).collect();
        match Frame::from_bytes(&frame_bytes) {
            Ok(frame) => Ok(Some(frame)),
            Err(e) => {
                log::warn!("Trama corrupta descartada: {}", e);
                Ok(None) // descartar trama corrupta, seguir intentando
            }
        }
    }
}
