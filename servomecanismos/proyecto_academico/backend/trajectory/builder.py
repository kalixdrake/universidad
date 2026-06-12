"""
Generación del trébol parametrizado, completa: trayectoria + aproximación Bezier.
"""

import numpy as np
from .trefoil import (
    generate_trefoil_raw,
    parametrize_constant_speed,
    find_tangential_entry,
    generate_bezier_approach,
    filter_savgol,
)


def build_full_trajectory(
    n_petals: int = 7,
    scale: float = 1.25,
    v_const_cm_s: float = 10.0,
    beta_deg: float = 45.0,
    L: float = 0.2,
    x_center: float = 0.2,
    y_center: float = 0.2,
    th1_home_deg: float = 120.0,
    th2_home_deg: float = -40.0,
    l1: float = 0.195,
    l2: float = 0.25,
    dt: float = 0.01,
) -> dict:
    """
    Construye la trayectoria completa del trébol con fase de aproximación.

    Retorna un dict con todas las señales calculadas.
    """
    # Convertir velocidad de cm/s a m/s
    v_const = v_const_cm_s / 100.0

    # ── 1. Trayectoria cruda ───────────────────────
    phi_raw, x_raw, y_raw, r_raw = generate_trefoil_raw(
        L=L, n_petals=n_petals, scale=scale, beta_deg=beta_deg
    )

    # ── 2. Parametrización a velocidad constante ───
    t, x_traj, y_traj, phi_t = parametrize_constant_speed(
        x_raw, y_raw, phi_raw,
        v_const=v_const, dt=dt,
        n_petals=n_petals, scale=scale, beta_deg=beta_deg,
        L=L, x_center=x_center, y_center=y_center,
    )

    # ── 3. Posición Home ───────────────────────────
    th1_home = np.deg2rad(th1_home_deg)
    th2_home = np.deg2rad(th2_home_deg)
    x_home = l1 * np.cos(th1_home) + l2 * np.cos(th1_home + th2_home)
    y_home = l1 * np.sin(th1_home) + l2 * np.sin(th1_home + th2_home)

    # ── 4. Entrada tangencial ──────────────────────
    x_traj_opt, y_traj_opt, idx_tang = find_tangential_entry(
        x_traj, y_traj, x_home, y_home
    )

    # ── 5. Curva de Bezier ─────────────────────────
    x_aprox, y_aprox = generate_bezier_approach(
        x_home, y_home,
        float(x_traj_opt[0]), float(y_traj_opt[0]),
        float(x_traj_opt[1]), float(y_traj_opt[1]),
        v_const=v_const, dt=dt,
    )

    # ── 6. Concatenar trayectoria completa ─────────
    x_total = np.concatenate([x_aprox, x_traj_opt])
    y_total = np.concatenate([y_aprox, y_traj_opt])

    return {
        'x_total': x_total,
        'y_total': y_total,
        'x_aprox': x_aprox,
        'y_aprox': y_aprox,
        'x_trebol': x_traj_opt,
        'y_trebol': y_traj_opt,
        'x_home': x_home,
        'y_home': y_home,
        'th1_home': th1_home,
        'th2_home': th2_home,
        'idx_tang': idx_tang,
        't_raw': t,
        'x_raw': x_traj,
        'y_raw': y_traj,
        'dt': dt,
        'L': L,
        'scale': scale,
        'x_center': x_center,
        'y_center': y_center,
    }
