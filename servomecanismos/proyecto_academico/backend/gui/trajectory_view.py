"""
Vista de monitoreo en tiempo real.

Muestra:
- Trayectoria de referencia (calculada) en azul
- Trayectoria real (telemetría ESP32) en rojo, dibujándose en vivo
- Gráficas de error de seguimiento
"""
import customtkinter as ctk
import numpy as np
import time
import threading
from collections import deque
from typing import Callable

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.lines import Line2D

from ..communication import TelemetryFrame


class TrajectoryView(ctk.CTkFrame):
    """Vista de monitoreo con gráficos en tiempo real."""

    MAX_TELEMETRY_PTS = 5000

    def __init__(self, parent, on_send: Callable, on_stop: Callable):
        super().__init__(parent)
        self._on_send = on_send
        self._on_stop = on_stop

        # Datos de telemetría
        self._telem_buf: deque[TelemetryFrame] = deque(maxlen=self.MAX_TELEMETRY_PTS)
        self._ref_trajectory: dict | None = None
        self._running = False
        self._start_time: float | None = None

        # ── Layout ────────────────────────────────────
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=0)   # botones
        self.grid_rowconfigure(1, weight=1)   # gráficos

        # Barra de botones
        self.btn_frame = ctk.CTkFrame(self, fg_color="transparent")
        self.btn_frame.grid(row=0, column=0, sticky="ew", pady=5, padx=5)

        self.btn_send = ctk.CTkButton(
            self.btn_frame, text="▶ Enviar a ESP32",
            font=ctk.CTkFont(weight="bold"),
            fg_color="#1A6B3C", hover_color="#145A30",
            command=self._send,
            state="disabled",
        )
        self.btn_send.grid(row=0, column=0, padx=5)

        self.btn_stop = ctk.CTkButton(
            self.btn_frame, text="⏹ Parada de Emergencia",
            fg_color="#8B1A1A", hover_color="#6B1414",
            command=self._stop,
        )
        self.btn_stop.grid(row=0, column=1, padx=5)

        self.lbl_status = ctk.CTkLabel(
            self.btn_frame, text="Esperando trayectoria...",
            font=ctk.CTkFont(size=12),
        )
        self.lbl_status.grid(row=0, column=2, padx=20, sticky="e")
        self.btn_frame.grid_columnconfigure(2, weight=1)

        # ── Figura Matplotlib ─────────────────────────
        self.fig = Figure(figsize=(10, 7), dpi=100, facecolor='#2B2B2B')

        # Subplots
        gs = self.fig.add_gridspec(2, 3, height_ratios=[3, 2])

        self.ax_cart = self.fig.add_subplot(gs[0, :2], facecolor='#2B2B2B')
        self.ax_error = self.fig.add_subplot(gs[1, :2], facecolor='#2B2B2B')
        self.ax_info = self.fig.add_subplot(gs[0, 2], facecolor='#2B2B2B')
        self.ax_state = self.fig.add_subplot(gs[1, 2], facecolor='#2B2B2B')

        self.fig.tight_layout(pad=2.5)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().grid(row=1, column=0, sticky="nsew", padx=5, pady=5)

        # Líneas (se crean en set_reference_trajectory)
        self._line_ref: Line2D | None = None
        self._line_real: Line2D | None = None
        self._line_error1: Line2D | None = None
        self._line_error2: Line2D | None = None

        # Timer de actualización
        self._update_interval = 50  # ms
        self._update_job: str | None = None

    # ── API ───────────────────────────────────────────
    def set_reference_trajectory(self, trajectory: dict):
        """Establece la trayectoria de referencia para dibujar."""
        self._ref_trajectory = trajectory
        self._telem_buf.clear()
        self.btn_send.configure(state="normal")
        self.lbl_status.configure(text="Listo para enviar")

        # Dibujar referencia
        self._draw_reference()
        self._start_animation()

    def on_telemetry(self, telem: TelemetryFrame):
        """Callback que recibe telemetría de la ESP32."""
        if self._start_time is None:
            self._start_time = time.time()

        self._telem_buf.append(telem)

        # Actualizar label
        state_names = {0: 'Idle', 1: 'Ejecutando', 2: 'Error', 3: 'Completado'}
        state_name = state_names.get(telem.state, '?')
        self.lbl_status.configure(
            text=f"Punto {telem.waypoint_index} | θ₁={np.rad2deg(telem.theta1):.1f}° "
                 f"θ₂={np.rad2deg(telem.theta2):.1f}° | {state_name}"
        )

        if telem.state == 3:
            self.lbl_status.configure(text="✓ Trayectoria completada")

    # ── Dibujo ────────────────────────────────────────
    def _draw_reference(self):
        """Dibuja la trayectoria de referencia."""
        if not self._ref_trajectory:
            return

        traj = self._ref_trajectory
        x_ref = traj.get('x_total', [])
        y_ref = traj.get('y_total', [])

        self.ax_cart.clear()
        self.ax_cart.plot(x_ref, y_ref, 'cyan', linewidth=0.6, alpha=0.5, label='Referencia')
        (self._line_ref,) = self.ax_cart.plot([], [], 'cyan', linewidth=1.2, label='Ref (progreso)')
        (self._line_real,) = self.ax_cart.plot([], [], 'r-', linewidth=1.5, label='Real (ESP32)')
        self.ax_cart.plot(x_ref[0], y_ref[0], 'go', markersize=5)
        self.ax_cart.set_title('Trayectoria Cartesiana', color='white')
        self.ax_cart.set_aspect('equal')
        self.ax_cart.legend(fontsize=7, loc='upper right')
        self.ax_cart.tick_params(colors='white', labelsize=7)
        for spine in self.ax_cart.spines.values():
            spine.set_color('#555555')

        self.ax_error.clear()
        self.ax_error.set_title('Error de Seguimiento [°]', color='white')
        self.ax_error.set_xlabel('Tiempo [s]', color='white')
        (self._line_error1,) = self.ax_error.plot([], [], 'orange', linewidth=0.8, label='e₁')
        (self._line_error2,) = self.ax_error.plot([], [], 'magenta', linewidth=0.8, label='e₂')
        self.ax_error.legend(fontsize=7)
        self.ax_error.tick_params(colors='white', labelsize=7)
        for spine in self.ax_error.spines.values():
            spine.set_color('#555555')

        self.ax_info.clear()
        self.ax_info.axis('off')
        self.ax_info.text(0.5, 0.5, 'Esperando telemetría...',
                          ha='center', va='center', color='#888888',
                          transform=self.ax_info.transAxes, fontsize=10)

        self.ax_state.clear()
        self.ax_state.axis('off')

        self.canvas.draw()

    def _update_plot(self):
        """Actualiza los gráficos con los datos de telemetría."""
        if not self._telem_buf or not self._ref_trajectory:
            self._update_job = self.after(self._update_interval, self._update_plot)
            return

        telem_list = list(self._telem_buf)
        n = len(telem_list)

        # Posiciones reales (cinemática directa)
        l1, l2 = 0.195, 0.25
        x_real = np.array([
            l1 * np.cos(t.theta1) + l2 * np.cos(t.theta1 + t.theta2)
            for t in telem_list
        ])
        y_real = np.array([
            l1 * np.sin(t.theta1) + l2 * np.sin(t.theta1 + t.theta2)
            for t in telem_list
        ])

        # Referencia (interpolar al mismo número de puntos)
        x_ref = np.array(self._ref_trajectory.get('x_total', []))
        y_ref = np.array(self._ref_trajectory.get('y_total', []))
        dyn = self._ref_trajectory.get('dynamics', {})
        theta1_ref = np.array(dyn.get('theta1', []))
        theta2_ref = np.array(dyn.get('theta2', []))

        # Actualizar trayectoria cartesiana
        if self._line_real:
            self._line_real.set_data(x_real, y_real)

        # Actualizar progreso de referencia
        wp = telem_list[-1].waypoint_index
        if self._line_ref and wp < len(x_ref):
            self._line_ref.set_data(x_ref[:wp+1], y_ref[:wp+1])
            self.ax_cart.relim()
            self.ax_cart.autoscale_view()

        # Error de seguimiento
        if self._line_error1 and len(theta1_ref) > 0 and n > 0:
            idx = min(wp, len(theta1_ref) - 1)
            e1 = np.rad2deg(np.array([t.theta1 for t in telem_list]) - theta1_ref[
                np.minimum(np.arange(n) + max(0, wp - n + 1), len(theta1_ref) - 1)
            ])
            e2 = np.rad2deg(np.array([t.theta2 for t in telem_list]) - theta2_ref[
                np.minimum(np.arange(n) + max(0, wp - n + 1), len(theta2_ref) - 1)
            ])
            t_arr = np.arange(n) * self._ref_trajectory.get('dt', 0.01) * 0.1
            self._line_error1.set_data(t_arr, e1)
            self._line_error2.set_data(t_arr, e2)
            self.ax_error.relim()
            self.ax_error.autoscale_view()

        # Info text
        if telem_list:
            last = telem_list[-1]
            self.ax_info.clear()
            self.ax_info.axis('off')
            info_lines = [
                f"Waypoint: {last.waypoint_index}",
                f"θ₁: {np.rad2deg(last.theta1):.1f}°",
                f"θ₂: {np.rad2deg(last.theta2):.1f}°",
                f"ω₁: {np.rad2deg(last.omega1):.1f}°/s",
                f"ω₂: {np.rad2deg(last.omega2):.1f}°/s",
                f"Estado: {['Idle','Run','Err','Done'][min(last.state,3)]}",
            ]
            self.ax_info.text(0.5, 0.5, '\n'.join(info_lines),
                              ha='center', va='center', color='white',
                              transform=self.ax_info.transAxes, fontsize=9,
                              fontfamily='monospace')

        self.canvas.draw_idle()
        self._update_job = self.after(self._update_interval, self._update_plot)

    def _start_animation(self):
        """Inicia el timer de actualización."""
        if self._update_job:
            self.after_cancel(self._update_job)
        self._update_job = self.after(self._update_interval, self._update_plot)

    # ── Botones ───────────────────────────────────────
    def _send(self):
        """Callback del botón Enviar."""
        self._telem_buf.clear()
        self._start_time = None
        self.btn_send.configure(state="disabled")
        self.lbl_status.configure(text="Enviando...")
        self._on_send()
        self.btn_send.configure(state="normal")

    def _stop(self):
        """Callback del botón Parada de Emergencia."""
        self._on_stop()
        self.lbl_status.configure(text="⚠ PARADA DE EMERGENCIA")
