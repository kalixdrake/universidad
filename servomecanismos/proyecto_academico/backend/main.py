#!/usr/bin/env python3
"""
Robot 2R — Backend de generación de trayectorias y comunicación con ESP32.

Uso:
    python -m backend.main

La aplicación:
1. Abre interfaz gráfica (customtkinter)
2. Busca/espera conexión con la ESP32 (red → USB)
3. Permite configurar parámetros del trébol con preview en tiempo real
4. Calcula cinemática inversa y dinámica
5. Envía trayectoria a la ESP32
6. Muestra telemetría en tiempo real
"""
import sys
import logging
from pathlib import Path

# Asegurar que backend/ está en el path
sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.gui import MainWindow

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    datefmt='%H:%M:%S',
)
logger = logging.getLogger(__name__)


def main():
    logger.info("Iniciando Robot 2R Backend...")
    app = MainWindow()
    app.mainloop()


if __name__ == "__main__":
    main()
