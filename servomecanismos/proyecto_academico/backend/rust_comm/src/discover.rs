//! Descubrimiento automático de ESP32 por USB y red.

use std::time::Duration;

/// Escanea puertos serie buscando una ESP32 por VID/PID o respuesta a comando ping.
pub fn discover_esp32_usb() -> Option<String> {
    let ports = serialport::available_ports().ok()?;

    for port_info in &ports {
        let port_name = &port_info.port_name;

        // Intento rápido: abrir, enviar ping, leer ACK
        if let Ok(mut port) = serialport::new(port_name, 115200)
            .timeout(Duration::from_millis(300))
            .open()
        {
            use std::io::{Read, Write};

            // Ping simple: enviamos bytes de sync y esperamos respuesta
            let ping = [0xAA, 0x55, 0x00, 0x00, 0x40]; // heartbeat type
            let _ = port.write_all(&ping);
            port.flush().ok();

            let mut buf = [0u8; 32];
            if let Ok(n) = port.read(&mut buf) {
                if n >= 7 && buf[0] == 0xAA && buf[1] == 0x55 {
                    log::info!("ESP32 detectada en {}", port_name);
                    return Some(port_name.clone());
                }
            }
        }
    }
    None
}

/// Busca la ESP32 en la red local vía UDP broadcast en puerto 4210.
pub fn discover_esp32_network() -> Option<String> {
    use socket2::{Domain, Protocol, Socket, Type};
    use std::net::{Ipv4Addr, SocketAddrV4, UdpSocket};

    let broadcast_addr = SocketAddrV4::new(Ipv4Addr::new(255, 255, 255, 255), 4210);
    let bind_addr = SocketAddrV4::new(Ipv4Addr::UNSPECIFIED, 0);

    let socket = UdpSocket::bind(bind_addr).ok()?;
    socket
        .set_broadcast(true)
        .ok()?;
    socket
        .set_read_timeout(Some(Duration::from_secs(2)))
        .ok()?;

    let discovery_msg = b"ESP32_DISCOVER";
    socket.send_to(discovery_msg, broadcast_addr).ok()?;

    let mut buf = [0u8; 128];
    match socket.recv_from(&mut buf) {
        Ok((n, src)) => {
            let response = String::from_utf8_lossy(&buf[..n]);
            if response.starts_with("ESP32_2R_ROBOT") {
                log::info!("ESP32 encontrada en red: {}", src.ip());
                // La respuesta incluye el puerto TCP: "ESP32_2R_ROBOT:8080"
                if let Some(port_str) = response.split(':').nth(1) {
                    if let Ok(port) = port_str.trim().parse::<u16>() {
                        return Some(format!("{}:{}", src.ip(), port));
                    }
                }
                Some(format!("{}:8080", src.ip()))
            } else {
                None
            }
        }
        Err(_) => None,
    }
}
