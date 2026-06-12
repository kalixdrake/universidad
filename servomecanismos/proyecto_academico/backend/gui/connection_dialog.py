"""
Diálogo de conexión con la ESP32.

Flujo:
1. Intenta conexión por red (WiFi)
2. Si falla, pide conectar por USB
3. Detecta ESP32, configura WiFi si es necesario
"""
import customtkinter as ctk
import threading
import logging
from typing import Callable

from ..communication.manager import CommunicationManager, ConnectionState

logger = logging.getLogger(__name__)


class ConnectionDialog(ctk.CTkFrame):
    """Panel de conexión con flujo guiado."""

    def __init__(self, parent, on_connected: Callable):
        super().__init__(parent)
        self._on_connected = on_connected
        self._comm_mgr = CommunicationManager()
        self._scanning = False

        # ── UI ────────────────────────────────────────
        self.grid_columnconfigure(0, weight=1)

        # Título
        self.title_label = ctk.CTkLabel(
            self, text="Conexión con ESP32",
            font=ctk.CTkFont(size=20, weight="bold"),
        )
        self.title_label.grid(row=0, column=0, pady=(20, 10))

        # Estado
        self.status_text = ctk.CTkTextbox(self, height=120, wrap="word")
        self.status_text.grid(row=1, column=0, sticky="ew", padx=20, pady=10)
        self.status_text.insert("1.0", "Buscando ESP32 en la red WiFi...\n")
        self.status_text.configure(state="disabled")

        # Barra de progreso
        self.progress = ctk.CTkProgressBar(self)
        self.progress.grid(row=2, column=0, sticky="ew", padx=20, pady=10)
        self.progress.set(0)
        self.progress.configure(mode="indeterminate")

        # Botones
        self.btn_frame = ctk.CTkFrame(self, fg_color="transparent")
        self.btn_frame.grid(row=3, column=0, pady=20)

        self.btn_scan = ctk.CTkButton(
            self.btn_frame, text="Reintentar búsqueda",
            command=self._start_scan,
        )
        self.btn_scan.grid(row=0, column=0, padx=10)

        self.btn_usb = ctk.CTkButton(
            self.btn_frame, text="Conectar por USB",
            command=self._connect_usb_manual,
            fg_color="#2B5B84",
        )
        self.btn_usb.grid(row=0, column=1, padx=10)

        # Manual USB entry
        self.usb_frame = ctk.CTkFrame(self, fg_color="transparent")
        self.usb_entry = ctk.CTkEntry(
            self.usb_frame, placeholder_text="/dev/ttyUSB0",
            width=200,
        )
        self.usb_entry.grid(row=0, column=0, padx=5)
        self.usb_btn = ctk.CTkButton(
            self.usb_frame, text="Conectar", width=80,
            command=self._connect_usb_entry,
        )
        self.usb_btn.grid(row=0, column=1, padx=5)

        # WiFi config (aparece después de USB)
        self.wifi_frame = ctk.CTkFrame(self, fg_color="transparent")
        self.wifi_label = ctk.CTkLabel(
            self.wifi_frame, text="Configurar WiFi en ESP32:",
            font=ctk.CTkFont(weight="bold"),
        )
        self.wifi_label.grid(row=0, column=0, columnspan=2, pady=(10, 5), sticky="w")

        self.ssid_entry = ctk.CTkEntry(
            self.wifi_frame, placeholder_text="SSID (nombre de red)",
            width=200,
        )
        self.ssid_entry.grid(row=1, column=0, padx=5, pady=5)

        self.pass_entry = ctk.CTkEntry(
            self.wifi_frame, placeholder_text="Contraseña", show="•",
            width=200,
        )
        self.pass_entry.grid(row=1, column=1, padx=5, pady=5)

        self.ip_entry = ctk.CTkEntry(
            self.wifi_frame, placeholder_text="IP estática (ej: 192.168.1.100)",
            width=200,
        )
        self.ip_entry.grid(row=2, column=0, padx=5, pady=5)
        self.ip_entry.insert(0, "192.168.4.1")

        self.wifi_btn = ctk.CTkButton(
            self.wifi_frame, text="Enviar configuración WiFi",
            command=self._configure_wifi,
        )
        self.wifi_btn.grid(row=2, column=1, padx=5, pady=5)

        # Iniciar escaneo automático
        self.after(500, self._start_scan)

    # ── Lógica ────────────────────────────────────────
    def _append_status(self, text: str):
        """Agrega texto al área de estado."""
        self.status_text.configure(state="normal")
        self.status_text.insert("end", text + "\n")
        self.status_text.see("end")
        self.status_text.configure(state="disabled")

    def _start_scan(self):
        """Inicia escaneo en hilo separado."""
        if self._scanning:
            return
        self._scanning = True
        self.progress.start()
        self._append_status("Escaneando...")
        threading.Thread(target=self._scan_thread, daemon=True).start()

    def _scan_thread(self):
        """Hilo de escaneo: red → USB."""
        # 1. Intentar red
        self.after(0, lambda: self._append_status("→ Buscando ESP32 en red WiFi (UDP broadcast)..."))
        net = self._comm_mgr.discover_network()
        if net:
            self.after(0, lambda: self._append_status(f"✓ ESP32 encontrada en {net}"))
            ok = self._comm_mgr.connect_network(net)
            if ok:
                self.after(0, lambda: self._on_connection_success())
                self.after(0, self.progress.stop)
                self._scanning = False
                return
            else:
                self.after(0, lambda: self._append_status("✗ Error conectando por red"))

        # 2. Intentar USB
        self.after(0, lambda: self._append_status("→ Buscando ESP32 en puertos USB..."))
        usb = self._comm_mgr.discover_usb()
        if usb:
            self.after(0, lambda: self._append_status(f"✓ ESP32 detectada en {usb}"))
            ok = self._comm_mgr.connect_usb(usb)
            if ok:
                self.after(0, lambda: self._append_status("✓ Conectado por USB. Puede configurar WiFi arriba."))
                self.after(0, self._show_wifi_config)
                self.after(0, lambda: self._on_connection_success())
                self.after(0, self.progress.stop)
                self._scanning = False
                return
            else:
                self.after(0, lambda: self._append_status("✗ Error conectando por USB"))

        # 3. No se encontró
        self.after(0, lambda: self._append_status(
            "⚠ No se encontró la ESP32.\n"
            "   Conecte la ESP32 por USB y presione 'Conectar por USB'\n"
            "   o ingrese el puerto manualmente."
        ))
        self.after(0, self.progress.stop)
        self._scanning = False

    def _connect_usb_manual(self):
        """Escanea puertos USB (sin pasar por red primero)."""
        self._append_status("Buscando en puertos USB...")
        threading.Thread(target=self._usb_thread, daemon=True).start()

    def _usb_thread(self):
        usb = self._comm_mgr.discover_usb()
        if usb:
            ok = self._comm_mgr.connect_usb(usb)
            if ok:
                self.after(0, lambda: self._append_status(f"✓ Conectado en {usb}"))
                self.after(0, self._show_wifi_config)
                self.after(0, lambda: self._on_connection_success())
                return
        self.after(0, lambda: self._append_status("✗ No se detectó ESP32 en USB"))

    def _connect_usb_entry(self):
        """Conecta al puerto ingresado manualmente."""
        port = self.usb_entry.get().strip()
        if not port:
            self._append_status("Ingrese un puerto válido")
            return
        ok = self._comm_mgr.connect_usb(port)
        if ok:
            self._append_status(f"✓ Conectado en {port}")
            self._show_wifi_config()
            self._on_connection_success()
        else:
            self._append_status(f"✗ No se pudo conectar a {port}")

    def _configure_wifi(self):
        """Envía configuración WiFi a la ESP32 por USB."""
        ssid = self.ssid_entry.get().strip()
        password = self.pass_entry.get().strip()
        static_ip = self.ip_entry.get().strip() or "192.168.4.1"

        if not ssid:
            self._append_status("Ingrese SSID")
            return

        ok = self._comm_mgr.configure_network(ssid, password, static_ip)
        if ok:
            self._append_status(f"✓ Configuración WiFi enviada: {ssid}")
        else:
            self._append_status("✗ Error enviando configuración WiFi")

    def _show_wifi_config(self):
        """Muestra el panel de configuración WiFi."""
        self.wifi_frame.grid(row=4, column=0, pady=10, padx=20, sticky="ew")
        self.usb_frame.grid(row=5, column=0, pady=5, padx=20, sticky="ew")

    def _on_connection_success(self):
        """Notifica a la ventana principal que la conexión fue exitosa."""
        self._on_connected(self._comm_mgr)
