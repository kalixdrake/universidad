# Robot 2R — Backend de Generación de Trayectorias

Backend en Python + Rust para la generación, cálculo cinemático/dinámico y comunicación en tiempo real con una ESP32 que controla un
manipulador robótico planar de 2 grados de libertad (2R).


## Arquitectura

```
┌──────────────────────────────────────────────────────────────────┐
│                        PC (Python + Rust)                         │
│  ┌────────────┐   ┌──────────────────┐   ┌────────────────────┐  │
│  │  GUI       │   │  Trayectorias    │   │  Comunicación      │  │
│  │ customtkinter│  │  numpy / scipy   │   │  pyserial / rust   │  │
│  │  • Sliders  │   │  • Trébol N-pétalos│  │  • USB (serial)   │  │
│  │  • Preview  │   │  • Bézier approach│  │  • WiFi (TCP/UDP)  │  │
│  │  • Monitoreo │  │  • Cinemática inv │   │  • Protocolo CRC32 │  │
│  └────────────┘   └──────────────────┘   └────────┬───────────┘  │
└───────────────────────────────────────────────────┬──────────────┘
                                                    │ USB / WiFi
┌───────────────────────────────────────────────────┴──────────────┐
│                        ESP32 (Rust firmware)                      │
│  ┌──────────────────┐   ┌──────────────┐   ┌──────────────────┐  │
│  │  Recepción       │   │  Control PID │   │  Telemetría      │  │
│  │  de trayectoria  │   │  2 motores DC│   │  100 Hz          │  │
│  └──────────────────┘   └──────────────┘   └──────────────────┘  │
└──────────────────────────────────────────────────────────────────┘
```

## Requisitos

### Sistema operativo
- Linux (Ubuntu/Debian), Windows 10/11, o macOS

---

### 🐧 Linux (Ubuntu/Debian)

#### Dependencias del sistema
```bash
sudo apt-get install -y python3 python3-venv python3-full \
    libudev-dev pkg-config build-essential curl
```

#### Instalar Rust
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source "$HOME/.cargo/env"
```

---

### 🪟 Windows 10/11

#### 1. Instalar Python

Descargar el instalador desde [python.org](https://www.python.org/downloads/) (versión 3.10 o superior).

> ⚠️ **Importante:** Durante la instalación, marcar la casilla **"Add Python to PATH"**.

Verificar la instalación:
```powershell
python --version
pip --version
```

#### 2. Instalar Rust

Descargar y ejecutar el instalador desde [rustup.rs](https://rustup.rs/).

O desde PowerShell:
```powershell
winget install Rustlang.Rustup
```

Verificar:
```powershell
rustc --version
cargo --version
```

#### 3. Opcional: Visual Studio Build Tools

Si al compilar `maturin develop` aparecen errores de linker, instalar
[Microsoft Visual C++ Build Tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/).

En el instalador, seleccionar:
- **"Desarrollo para el escritorio con C++"** (Desktop development with C++)
- Asegurar que incluya **Windows SDK** y **MSVC v143**

> **Nota:** En la mayoría de los casos, Rust ya incluye su propio linker y no
> necesitas instalar Visual Studio. Solo es necesario si `maturin` falla con
> errores tipo `link.exe not found`.

---

## Instalación y ejecución

### 1. Clonar el repositorio

```bash
git clone <url-del-repo>
cd proyecto_academico
```

### 2. Crear entorno virtual e instalar dependencias Python

**Linux / macOS:**
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r backend/requirements.txt
```

**Windows (PowerShell / CMD):**
```powershell
python -m venv .venv
.venv\Scripts\activate
pip install -r backend\requirements.txt
```

### 3. Compilar e instalar el módulo Rust (opcional pero recomendado)

```bash
cd backend/rust_comm

# Compilar como extensión nativa de Python
maturin develop --release

cd ../..
```

> **Nota:** Si no compilas el módulo Rust, la aplicación usará automáticamente
> `pyserial` + `socket` como fallback. La funcionalidad es idéntica, solo
> cambia el rendimiento de la comunicación serie.

### 4. Ejecutar la aplicación

**Linux / macOS:**
```bash
.venv/bin/python backend/main.py
```

**Windows:**
```powershell
.venv\Scripts\python backend\main.py
```

## Uso de la interfaz

| Pestaña | Función |
|---------|---------|
| **1. Conexión** | Busca la ESP32 en red WiFi o USB. Permite configurar WiFi en la ESP32 vía USB. |
| **2. Configuración** | Sliders interactivos con preview en tiempo real del trébol, torques y ángulos. |
| **3. Monitoreo** | Vista de la trayectoria de referencia (azul) y la trayectoria real (roja) recibida por telemetría. |

### Parámetros del trébol

| Parámetro | Rango | Descripción |
|-----------|-------|-------------|
| Pétalos | 3 – 7 | Número de hojas del trébol |
| Escala | 0.75 – 1.25 | Factor de escala de la figura |
| Velocidad | 1 – 10 cm/s | Velocidad tangencial de dibujado |
| Rotación | ±45° | Ángulo de rotación del trébol |
| Ciclos | 1 – 10 | Veces que se repite la trayectoria |

## Estructura del proyecto

```
proyecto_academico/
├── backend/
│   ├── main.py                      # Punto de entrada
│   ├── requirements.txt             # Dependencias Python
│   │
│   ├── gui/                         # Interfaz gráfica (customtkinter)
│   │   ├── main_window.py           # Ventana principal con tabs
│   │   ├── connection_dialog.py     # Flujo conexión red → USB
│   │   ├── config_panel.py          # Sliders + preview interactivo
│   │   └── trajectory_view.py       # Monitoreo en tiempo real
│   │
│   ├── trajectory/                  # Cálculo de trayectorias
│   │   ├── trefoil.py               # Generación del trébol + Bézier
│   │   ├── builder.py               # Orquestador de trayectoria
│   │   └── kinematics.py            # Cinemática inversa + dinámica
│   │
│   ├── communication/               # Protocolo PC ↔ ESP32
│   │   ├── __init__.py              # Trama binaria + CRC32
│   │   └── manager.py               # Gestor dual (Rust / pyserial)
│   │
│   ├── rust_comm/                   # Módulo Rust (PyO3) — PC
│   │   ├── Cargo.toml
│   │   ├── README.md
│   │   └── src/
│   │       ├── lib.rs               # Esp32Connection (Python API)
│   │       ├── protocol.rs          # Protocolo binario con CRC32
│   │       ├── discover.rs          # Auto-detección USB + UDP
│   │       └── serial_link.rs       # Enlace serie con buffering
│   │
│   └── esp32_firmware/              # Firmware ESP32 (Rust no_std)
│       ├── Cargo.toml
│       ├── README.md
│       └── src/
│           └── main.rs              # Control PID + comunicaciones
│
├── caracterizacion_trayectoria_trebol.m   # MATLAB original
├── modelo_motor_dc.m
├── caso_critico.m
├── datasheets/
└── esp32/                           # Código ESP32 (PlatformIO/C++)
```

## Licencia

Proyecto académico — Universidad, Servomecanismos.
