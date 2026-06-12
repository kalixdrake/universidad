//! Protocolo binario empaquetado con CRC32 para comunicación PC ↔ ESP32.
//!
//! Formato de trama:
//! ┌──────────┬──────────┬────────────┬───────────────────┬──────────┐
//! │ 0xAA 0x55│  length  │  msg_type  │  payload (N bytes) │  CRC32   │
//! │  2 bytes │  2 bytes │   1 byte   │     N bytes        │  4 bytes │
//! └──────────┴──────────┴────────────┴───────────────────┴──────────┘

use serde::{Deserialize, Serialize};

pub const SYNC_WORD: [u8; 2] = [0xAA, 0x55];
pub const HEADER_SIZE: usize = 5; // sync(2) + length(2) + type(1)
pub const CRC_SIZE: usize = 4;
pub const MAX_PAYLOAD: usize = 16 * 1024; // 16 KB máximo por trama

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Frame {
    /// Configuración de red WiFi
    NetworkConfig(NetworkConfig),
    /// Datos de trayectoria (puede fragmentarse en múltiples tramas)
    Trajectory(TrajectoryData),
    /// Segmento de trayectoria (para fragmentación)
    TrajectoryChunk(TrajectoryChunk),
    /// Comando de inicio de ejecución
    StartExecution,
    /// Comando de parada de emergencia
    EmergencyStop,
    /// ACK / NACK
    Ack { success: bool, message: String },
    /// Telemetría en tiempo real desde la ESP32
    Telemetry(TelemetryFrame),
    /// Heartbeat / keep-alive
    Heartbeat { seq: u32 },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NetworkConfig {
    pub ssid: String,
    pub password: String,
    pub static_ip: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrajectoryData {
    pub theta1: Vec<f64>,
    pub theta2: Vec<f64>,
    pub omega1: Vec<f64>,
    pub omega2: Vec<f64>,
    pub dt: f64,
    pub cycles: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrajectoryChunk {
    pub chunk_index: u32,
    pub total_chunks: u32,
    pub theta1: Vec<f64>,
    pub theta2: Vec<f64>,
    pub omega1: Vec<f64>,
    pub omega2: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TelemetryFrame {
    /// Timestamp en microsegundos desde el inicio
    pub timestamp_us: u64,
    /// Posición actual θ1 [rad]
    pub theta1: f64,
    /// Posición actual θ2 [rad]
    pub theta2: f64,
    /// Velocidad actual ω1 [rad/s]
    pub omega1: f64,
    /// Velocidad actual ω2 [rad/s]
    pub omega2: f64,
    /// Índice del punto de trayectoria actual
    pub waypoint_index: u32,
    /// Estado: 0=idle, 1=running, 2=error, 3=completed
    pub state: u8,
    /// Código de error (si state=2)
    pub error_code: u8,
}

impl Frame {
    /// Serializa la trama a bytes con CRC32.
    pub fn to_bytes(&self) -> Vec<u8> {
        let payload = bincode::serialize(self).unwrap_or_default();
        let payload_len = payload.len();
        let total_len = HEADER_SIZE + payload_len + CRC_SIZE;

        let mut buf = Vec::with_capacity(total_len);
        buf.extend_from_slice(&SYNC_WORD);
        buf.extend_from_slice(&(payload_len as u16).to_le_bytes());
        buf.extend_from_slice(&[self.msg_type()]);
        buf.extend_from_slice(&payload);

        let crc = crc32fast::hash(&buf[..HEADER_SIZE + payload_len]);
        buf.extend_from_slice(&crc.to_le_bytes());

        buf
    }

    /// Deserializa una trama desde bytes (debe incluir CRC).
    pub fn from_bytes(data: &[u8]) -> Result<Self, String> {
        if data.len() < HEADER_SIZE + CRC_SIZE {
            return Err("Trama demasiado corta".into());
        }
        if data[0] != SYNC_WORD[0] || data[1] != SYNC_WORD[1] {
            return Err("Sync word inválida".into());
        }

        let payload_len = u16::from_le_bytes([data[2], data[3]]) as usize;
        let msg_type = data[4];

        if data.len() < HEADER_SIZE + payload_len + CRC_SIZE {
            return Err("Longitud de trama insuficiente".into());
        }

        let payload_end = HEADER_SIZE + payload_len;
        let expected_crc = u32::from_le_bytes([
            data[payload_end],
            data[payload_end + 1],
            data[payload_end + 2],
            data[payload_end + 3],
        ]);
        let actual_crc = crc32fast::hash(&data[..payload_end]);

        if actual_crc != expected_crc {
            return Err(format!(
                "CRC mismatch: esperado 0x{:08X}, calculado 0x{:08X}",
                expected_crc, actual_crc
            ));
        }

        bincode::deserialize(&data[HEADER_SIZE..payload_end])
            .map_err(|e| format!("Error deserializando: {}", e))
    }

    fn msg_type(&self) -> u8 {
        match self {
            Frame::NetworkConfig(_) => 0x01,
            Frame::Trajectory(_) => 0x02,
            Frame::TrajectoryChunk(_) => 0x03,
            Frame::StartExecution => 0x10,
            Frame::EmergencyStop => 0x11,
            Frame::Ack { .. } => 0x20,
            Frame::Telemetry(_) => 0x30,
            Frame::Heartbeat { .. } => 0x40,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_roundtrip_telemetry() {
        let original = Frame::Telemetry(TelemetryFrame {
            timestamp_us: 12345678,
            theta1: 1.57,
            theta2: -0.78,
            omega1: 0.1,
            omega2: -0.05,
            waypoint_index: 42,
            state: 1,
            error_code: 0,
        });
        let bytes = original.to_bytes();
        let decoded = Frame::from_bytes(&bytes).unwrap();
        match (&original, &decoded) {
            (Frame::Telemetry(a), Frame::Telemetry(b)) => {
                assert!((a.theta1 - b.theta1).abs() < 1e-9);
                assert!((a.theta2 - b.theta2).abs() < 1e-9);
                assert_eq!(a.waypoint_index, b.waypoint_index);
            }
            _ => panic!("Tipo incorrecto"),
        }
    }

    #[test]
    fn test_crc_failure() {
        let frame = Frame::Ack {
            success: true,
            message: "ok".into(),
        };
        let mut bytes = frame.to_bytes();
        // Corromper el CRC
        let len = bytes.len();
        bytes[len - 1] ^= 0xFF;
        assert!(Frame::from_bytes(&bytes).is_err());
    }
}
