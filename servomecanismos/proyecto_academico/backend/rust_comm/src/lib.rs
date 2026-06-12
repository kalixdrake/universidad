//! rust_comm — Módulo de comunicación PC ↔ ESP32
//!
//! Proporciona descubrimiento automático, conexión serie/USB, y protocolo
//! empaquetado con verificación CRC32 para intercambio fiable de datos de
//! trayectoria entre el backend Python y la ESP32.

mod discover;
mod protocol;
mod serial_link;

use pyo3::prelude::*;
use std::sync::Mutex;

use discover::{discover_esp32_usb, discover_esp32_network};
use protocol::Frame;
use serial_link::SerialLink;

/// Datos de telemetría expuestos a Python.
#[pyclass]
#[derive(Clone, Debug)]
pub struct TelemetryData {
    #[pyo3(get)]
    pub timestamp_us: u64,
    #[pyo3(get)]
    pub theta1: f64,
    #[pyo3(get)]
    pub theta2: f64,
    #[pyo3(get)]
    pub omega1: f64,
    #[pyo3(get)]
    pub omega2: f64,
    #[pyo3(get)]
    pub waypoint_index: u32,
    #[pyo3(get)]
    pub state: u8,
    #[pyo3(get)]
    pub error_code: u8,
}

/// Estado de la conexión con la ESP32.
#[pyclass]
pub struct Esp32Connection {
    link: Mutex<Option<SerialLink>>,
    network_addr: Mutex<Option<String>>,
}

#[pymethods]
impl Esp32Connection {
    #[new]
    fn new() -> Self {
        Esp32Connection {
            link: Mutex::new(None),
            network_addr: Mutex::new(None),
        }
    }

    /// Escanea puertos USB en busca de una ESP32.
    /// Retorna el puerto (ej: "/dev/ttyUSB0") o None.
    fn discover_usb(&self) -> PyResult<Option<String>> {
        Ok(discover_esp32_usb())
    }

    /// Escanea la red local en busca de la ESP32 vía mDNS/UDP broadcast.
    /// Retorna dirección IP:puerto o None.
    fn discover_network(&self) -> PyResult<Option<String>> {
        Ok(discover_esp32_network())
    }

    /// Conecta a la ESP32 por USB (serial).
    fn connect_usb(&self, port: &str, baudrate: u32) -> PyResult<bool> {
        match SerialLink::open(port, baudrate) {
            Ok(link) => {
                *self.link.lock().unwrap() = Some(link);
                Ok(true)
            }
            Err(e) => Err(pyo3::exceptions::PyConnectionError::new_err(format!(
                "Error al conectar USB {}: {}",
                port, e
            ))),
        }
    }

    /// Envía configuración de red a la ESP32 por USB para habilitar WiFi.
    fn send_network_config(&self, ssid: &str, password: &str, static_ip: &str) -> PyResult<bool> {
        let cfg = protocol::NetworkConfig {
            ssid: ssid.to_string(),
            password: password.to_string(),
            static_ip: static_ip.to_string(),
        };
        let frame = Frame::NetworkConfig(cfg);
        self.send_frame(&frame)
    }

    /// Envía un vector de posiciones articulares (θ1, θ2) y velocidades.
    fn send_trajectory(
        &self,
        theta1: Vec<f64>,
        theta2: Vec<f64>,
        omega1: Vec<f64>,
        omega2: Vec<f64>,
        dt: f64,
        cycles: u32,
    ) -> PyResult<bool> {
        let traj = protocol::TrajectoryData {
            theta1,
            theta2,
            omega1,
            omega2,
            dt,
            cycles,
        };
        let frame = Frame::Trajectory(traj);
        self.send_frame(&frame)
    }

    /// Lee un frame de telemetría desde la ESP32 (bloqueante con timeout).
    fn read_telemetry(&self, timeout_ms: u64) -> PyResult<Option<TelemetryData>> {
        let mut link_opt = self.link.lock().unwrap();
        match link_opt.as_mut() {
            Some(link) => match link.read_frame(timeout_ms) {
                Ok(Some(frame)) => {
                    if let Frame::Telemetry(t) = frame {
                        Ok(Some(TelemetryData {
                            timestamp_us: t.timestamp_us,
                            theta1: t.theta1,
                            theta2: t.theta2,
                            omega1: t.omega1,
                            omega2: t.omega2,
                            waypoint_index: t.waypoint_index,
                            state: t.state,
                            error_code: t.error_code,
                        }))
                    } else {
                        Ok(None)
                    }
                }
                Ok(None) => Ok(None),
                Err(e) => Err(pyo3::exceptions::PyIOError::new_err(format!(
                    "Error lectura: {}",
                    e
                ))),
            },
            None => Err(pyo3::exceptions::PyConnectionError::new_err(
                "No hay conexión activa",
            )),
        }
    }

    /// Cierra la conexión activa.
    fn close(&self) {
        *self.link.lock().unwrap() = None;
    }

    /// Verifica si hay conexión activa.
    fn is_connected(&self) -> bool {
        self.link.lock().unwrap().is_some()
    }
}

// ── Internos ──────────────────────────────────────────────
impl Esp32Connection {
    fn send_frame(&self, frame: &Frame) -> PyResult<bool> {
        let mut link_opt = self.link.lock().unwrap();
        match link_opt.as_mut() {
            Some(link) => match link.write_frame(frame) {
                Ok(_) => Ok(true),
                Err(e) => Err(pyo3::exceptions::PyIOError::new_err(format!(
                    "Error envío: {}",
                    e
                ))),
            },
            None => Err(pyo3::exceptions::PyConnectionError::new_err(
                "No hay conexión activa",
            )),
        }
    }
}

// ── Módulo Python ─────────────────────────────────────────
#[pymodule]
fn rust_comm(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Esp32Connection>()?;
    m.add_class::<TelemetryData>()?;
    Ok(())
}
