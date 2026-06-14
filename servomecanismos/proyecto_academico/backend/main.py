#!/usr/bin/env python3
"""
Punto de entrada principal para el control y calibración del Robot 2R en la PC.
"""

import sys
from trajectory.generation import generate_clover_trajectory
from kinematics.inverse import inverse_kinematics
from communication.serial_conn import SerialConnection
from calibration.calibration_manager import CalibrationManager

def main():
    print("=" * 60)
    print("   SISTEMA DE CONTROL Y CALIBRACIÓN - ROBOT PLANAR 2R   ")
    print("=" * 60)
    
    # 1. Comprobación del Sistema (Esqueleto)
    print("\n[INFO] Inicializando módulos del PC...")
    
    # Ejemplo de llamada barebone a generación de trayectorias
    print("[INFO] Probando generador de trayectorias...")
    puntos = generate_clover_trajectory(petals=3, num_points=10)
    print(f"       -> Generados {len(puntos)} puntos de la trayectoria.")
    
    # Ejemplo de llamada barebone a cinemática inversa
    print("[INFO] Probando resolvedor de cinemática inversa...")
    theta1, theta2 = inverse_kinematics(x=0.1, y=0.2)
    print(f"       -> Ángulos calculados: theta1 = {theta1:.2f}°, theta2 = {theta2:.2f}°")
    
    # Ejemplo de inicialización de gestor de calibración
    print("[INFO] Cargando manejador de calibración...")
    manager = CalibrationManager()
    status = manager.get_status()
    print(f"       -> Estado actual de calibración: {status}")
    
    print("\n[INFO] Esqueleto del backend inicializado correctamente.")
    print("       Modifica los módulos bajo 'backend/' para añadir la funcionalidad real.")
    print("=" * 60)

if __name__ == "__main__":
    main()
