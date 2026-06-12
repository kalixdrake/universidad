# Módulo de Comunicación Rust (PyO3)

Este crate proporciona comunicación serie de alto rendimiento entre Python y la ESP32, compilado como extensión nativa de Python vía PyO3.

## Compilación

```bash
# Instalar maturin (solo una vez)
pip install maturin

# Compilar en modo release
cd backend/rust_comm
maturin develop --release
```

Esto compila el crate y lo instala como módulo `rust_comm` en tu entorno Python actual.

## Uso desde Python

```python
import rust_comm

conn = rust_comm.Esp32Connection()

# Descubrir ESP32
port = conn.discover_usb()
if port:
    conn.connect_usb(port, 115200)

# Enviar trayectoria
conn.send_trajectory(theta1, theta2, omega1, omega2, dt=0.01, cycles=1)

# Leer telemetría
telem = conn.read_telemetry(timeout_ms=100)
if telem:
    print(f"θ1={telem.theta1}, θ2={telem.theta2}")

conn.close()
```

## Estructura

```
src/
├── lib.rs         # Punto de entrada PyO3, clase Esp32Connection
├── protocol.rs    # Protocolo de tramas con CRC32
├── discover.rs    # Descubrimiento automático (USB + UDP)
└── serial_link.rs # Enlace serie con buffering
```

## Dependencias

- `pyo3` — Bindings Python ↔ Rust
- `serialport` — Comunicación serie nativa
- `socket2` — Sockets de red avanzados
- `bincode` + `serde` — Serialización binaria
- `crc32fast` — Integridad de tramas
