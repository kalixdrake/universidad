"""
Panel de configuración interactivo del trébol.

Muestra sliders para:
- Número de pétalos (3-7)
- Escala (0.75-1.25)
- Velocidad tangencial (1-10 cm/s)
- Rotación (±45°)
- Número de ciclos (1-10)

Con preview en tiempo real de la trayectoria y cinemática.
"""
import customtkinter as ctk
import numpy as np
import threading
import logging
from typing import Callable

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from ..trajectory import build_full_trajectory
from ..trajectory.kinematics import inverse_kinematics, inverse_dynamics, forward_kinematics

logger = logging.getLogger(__name__)


class ConfigPanel(ctk.CTkFrame):
    """Panel de configuración con preview interactivo."""

    def __init__(self, parent, on_calculate: Callable):
        super().__init__(parent)
        self._on_calculate = on_calculate
        self._calc_lock = threading.Lock()
        self._pending_calc: str | None = None

        # ── Layout: dos columnas ──────────────────────
        self.grid_columnconfigure(0, weight=0)  # controles
        self.grid_columnconfigure(1, weight=1)  # preview
        self.grid_rowconfigure(0, weight=1)

        # ── Columna izquierda: Controles ──────────────
        self.ctrl_frame = ctk.CTkScrollableFrame(self, width=320)
        self.ctrl_frame.grid(row=0, column=0, sticky="ns", padx=10, pady=10)

        self._build_controls()

        # ── Columna derecha: Preview ──────────────────
        self.preview_frame = ctk.CTkFrame(self)
        self.preview_frame.grid(row=0, column=1, sticky="nsew", padx=10, pady=10)

        # Matplotlib figure
        self.fig = Figure(figsize=(6, 5), dpi=100, facecolor='#2B2B2B')
        self.ax_traj = self.fig.add_subplot(2, 2, 1, facecolor='#2B2B2B')
        self.ax_torque = self.fig.add_subplot(2, 2, 2, facecolor='#2B2B2B')
        self.ax_angles = self.fig.add_subplot(2, 1, 2, facecolor='#2B2B2B')
        self.fig.tight_layout(pad=2.0)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.preview_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        # ── Botón Calcular ────────────────────────────
        self.btn_calc = ctk.CTkButton(
            self, text="▶ Calcular y enviar trayectoria",
            font=ctk.CTkFont(size=14, weight="bold"),
            height=40,
            command=self._calculate,
        )
        self.btn_calc.grid(row=1, column=0, columnspan=2, pady=10, padx=10, sticky="ew")

        # Hacer primer preview
        self.after(300, self._update_preview)

    def _build_controls(self):
        """Construye los sliders y labels."""
        row = 0

        # Título
        ctk.CTkLabel(
            self.ctrl_frame, text="Parámetros del Trébol",
            font=ctk.CTkFont(size=16, weight="bold"),
        ).grid(row=row, column=0, columnspan=2, pady=(0, 10), sticky="w")
        row += 1

        # Pétalos
        self._make_slider(row, "Número de pétalos", "petals", 3, 7, 7, is_int=True)
        row += 2

        # Escala
        self._make_slider(row, "Escala", "scale", 0.75, 1.25, 1.25, step=0.01)
        row += 2

        # Velocidad
        self._make_slider(row, "Velocidad tangencial (cm/s)", "speed", 1.0, 10.0, 10.0, step=0.1)
        row += 2

        # Rotación
        self._make_slider(row, "Rotación (°)", "rotation", -45.0, 45.0, 45.0, step=1.0)
        row += 2

        # Ciclos
        self._make_slider(row, "Número de ciclos", "cycles", 1, 10, 1, is_int=True)
        row += 3

        # Labels informativos
        self.info_label = ctk.CTkLabel(
            self.ctrl_frame, text="", justify="left",
            font=ctk.CTkFont(size=11),
        )
        self.info_label.grid(row=row, column=0, columnspan=2, pady=(20, 0), sticky="w")

    def _make_slider(
        self, row: int, label: str, key: str,
        min_val: float, max_val: float, default: float,
        step: float = 1.0, is_int: bool = False,
    ):
        """Crea un slider con label y value display."""
        var = ctk.DoubleVar(value=default)
        setattr(self, f'_var_{key}', var)

        ctk.CTkLabel(self.ctrl_frame, text=label).grid(
            row=row, column=0, columnspan=2, sticky="w", pady=(5, 0)
        )

        val_label = ctk.CTkLabel(self.ctrl_frame, text=str(default), width=40)
        val_label.grid(row=row + 1, column=1, sticky="e", padx=(5, 0))
        setattr(self, f'_lbl_{key}', val_label)

        slider = ctk.CTkSlider(
            self.ctrl_frame,
            from_=min_val, to=max_val,
            number_of_steps=int((max_val - min_val) / step),
            variable=var,
            command=lambda v, k=key, il=is_int: self._on_slider_change(k, v, il),
        )
        slider.grid(row=row + 1, column=0, sticky="ew", pady=(0, 5))
        setattr(self, f'_sld_{key}', slider)

    def _on_slider_change(self, key: str, value: float, is_int: bool):
        """Actualiza el label y agenda un preview."""
        lbl = getattr(self, f'_lbl_{key}')
        if is_int:
            lbl.configure(text=str(int(value)))
        else:
            lbl.configure(text=f"{value:.2f}")

        # Debounce: agenda preview después de 200ms
        if self._pending_calc:
            self.after_cancel(self._pending_calc)
        self._pending_calc = self.after(200, self._update_preview)

    def _get_params(self) -> dict:
        """Obtiene los parámetros actuales de los sliders."""
        return {
            'n_petals': int(self._var_petals.get()),
            'scale': self._var_scale.get(),
            'v_const_cm_s': self._var_speed.get(),
            'beta_deg': self._var_rotation.get(),
            'cycles': int(self._var_cycles.get()),
        }

    def _update_preview(self):
        """Calcula y dibuja el preview en tiempo real."""
        self._pending_calc = None
        params = self._get_params()

        def compute():
            try:
                traj = build_full_trajectory(
                    n_petals=params['n_petals'],
                    scale=params['scale'],
                    v_const_cm_s=params['v_const_cm_s'],
                    beta_deg=params['beta_deg'],
                )
                theta1_raw, theta2_raw = inverse_kinematics(
                    traj['x_total'], traj['y_total'],
                    traj['th1_home'], traj['th2_home'],
                )
                dyn = inverse_dynamics(theta1_raw, theta2_raw, traj['dt'])
                traj['dynamics'] = dyn
                traj['cycles'] = params['cycles']
                self.after(0, lambda: self._draw_preview(traj, dyn, params))
            except Exception as e:
                logger.error(f"Error en preview: {e}")

        threading.Thread(target=compute, daemon=True).start()

    def _draw_preview(self, traj: dict, dyn: dict, params: dict):
        """Dibuja los gráficos de preview."""
        self.ax_traj.clear()
        self.ax_torque.clear()
        self.ax_angles.clear()

        x_t, y_t = traj['x_total'], traj['y_total']

        # Trayectoria cartesiana
        self.ax_traj.plot(x_t, y_t, 'cyan', linewidth=0.8)
        self.ax_traj.plot(x_t[0], y_t[0], 'go', markersize=4, label='Inicio')
        self.ax_traj.plot(x_t[-1], y_t[-1], 'ro', markersize=4, label='Fin')
        self.ax_traj.set_title('Trayectoria Cartesiana', color='white', fontsize=9)
        self.ax_traj.set_aspect('equal')
        self.ax_traj.tick_params(colors='white', labelsize=7)
        self.ax_traj.legend(fontsize=6, loc='upper right')
        for spine in self.ax_traj.spines.values():
            spine.set_color('#555555')

        # Torques
        tau1 = dyn['tau1']
        tau2 = dyn['tau2']
        t_vec = np.arange(len(tau1)) * traj['dt']
        self.ax_torque.plot(t_vec, tau1, 'orange', linewidth=0.8, label='τ₁')
        self.ax_torque.plot(t_vec, tau2, 'magenta', linewidth=0.8, label='τ₂')
        self.ax_torque.set_title('Torques [N·m]', color='white', fontsize=9)
        self.ax_torque.tick_params(colors='white', labelsize=7)
        self.ax_torque.legend(fontsize=6)
        for spine in self.ax_torque.spines.values():
            spine.set_color('#555555')

        # Ángulos
        self.ax_angles.plot(t_vec, np.rad2deg(dyn['theta1']), 'cyan', linewidth=0.8, label='θ₁')
        self.ax_angles.plot(t_vec, np.rad2deg(dyn['theta2']), 'yellow', linewidth=0.8, label='θ₂')
        self.ax_angles.set_title('Ángulos Articulares [°]', color='white', fontsize=9)
        self.ax_angles.set_xlabel('Tiempo [s]', color='white', fontsize=8)
        self.ax_angles.tick_params(colors='white', labelsize=7)
        self.ax_angles.legend(fontsize=6)
        for spine in self.ax_angles.spines.values():
            spine.set_color('#555555')

        self.fig.tight_layout(pad=2.0)
        self.canvas.draw()

        # Actualizar info
        total_pts = len(traj['x_total'])
        duration = total_pts * traj['dt']
        max_tau = max(max(abs(tau1)), max(abs(tau2)))
        info = (
            f"Pétalos: {params['n_petals']} | Escala: {params['scale']:.2f} | "
            f"Vel: {params['v_const_cm_s']:.1f} cm/s | Rot: {params['beta_deg']:.0f}°\n"
            f"Puntos: {total_pts} | Duración: {duration:.1f}s | "
            f"Ciclos: {params['cycles']}\n"
            f"Torque máx: {max_tau:.3f} N·m"
        )
        self.info_label.configure(text=info)

        # Guardar trayectoria completa para envío
        self._cached_trajectory = traj

    def _calculate(self):
        """Callback del botón Calcular."""
        if not hasattr(self, '_cached_trajectory') or self._cached_trajectory is None:
            self._update_preview()
            self.after(500, self._calculate)
            return
        self._on_calculate(self._cached_trajectory)
