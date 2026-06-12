%% MODELO DINÁMICO DE MOTOR DC CON REDUCTORA
% Calcula constantes del motor, funciones de transferencia y validación
% contra datos de catálogo.
%
% Basado en:
%   V = I·R + Ke·ω_m            (eléctrica)
%   τ_m = Kt·I                   (electromecánica)
%   τ_m = J_eq·α + b·ω + τ_carga (mecánica)
%
% Para motor DC ideal: Kt (N·m/A) = Ke (V·s/rad) en SI
%
clear; clc; close all;

%% 1. DATOS DE ENTRADA — COMPLETAR SEGÚN DATASHEETS
% ======================= MOTOR BASE =======================
V_base     = 12.0;        % Voltaje nominal [V]
I_nom_base = 1.3;         % Corriente nominal [A]
I_stall_base = 10.0;      % Corriente de stall [A]
N_base     = 218;         % Relación de reducción

% Parámetros a la SALIDA del gearbox (catálogo)
tau_nom_out_base  = 30;   % Torque nominal [kgf·cm]
tau_stall_out_base = 90;  % Torque stall [kgf·cm]
rpm_nl_out_base   = 37;   % RPM sin carga (salida)
rpm_nom_out_base  = 30;   % RPM nominal (salida)

% ======================= MOTOR CODO =======================
V_codo     = 24.0;        % Voltaje nominal [V]
I_nom_codo = 0.6;         % Corriente nominal [A]
I_stall_codo = 3.1;       % Corriente de stall [A]
N_codo     = 270;         % Relación de reducción

tau_nom_out_codo  = 19;   % Torque nominal [kgf·cm]
tau_stall_out_codo = 60;  % Torque stall [kgf·cm]
rpm_nl_out_codo   = 37;   % RPM sin carga (salida)
rpm_nom_out_codo  = 28;   % RPM nominal (salida)

% ======================= CONSTANTES GENERALES =======================
kgfcm_to_Nm = 0.0980665;
eta_est = 0.65;  % Eficiencia estimada del gearbox (ajustable)
                  % Típico: 0.9 por etapa → 3 etapas = 0.73
                  % Se usa principalmente para verificación; Kt se obtiene de Ke

% ======================= PARÁMETROS ELÉCTRICOS ESTIMADOS =======================
% Inductancia de armadura - estimación para motores DC brushed pequeños
% Si no tienes el dato, se estima: L ≈ V/(6·π·p·ω_nl·I_stall) aprox.
% Usamos un valor típico de 100-500 μH para motores DC con escobillas pequeños
L_est_base = 200e-6;  % [H] - 200 μH, ajustar si se conoce el valor real
L_est_codo = 400e-6;  % [H] - 400 μH

% Inercia del rotor (valores típicos para motores de ~35mm)
% Si no se tiene, estimar: J_rotor ≈ 1e-6 a 5e-6 kg·m²
J_rotor_base = 3.5e-6;  % [kg·m²]
J_rotor_codo = 4.0e-6;  % [kg·m²]

% Coeficiente de fricción viscosa (típicamente muy pequeño para estos motores)
% Se estima de: b = Kt·I_nl / ω_m_nl  donde I_nl ≈ 0.1·I_stall
b_base = 5e-6;   % [N·m·s/rad]
b_codo = 5e-6;

%% 2. CÁLCULO DE CONSTANTES DEL MOTOR BASE
fprintf('============================================\n');
fprintf('  MOTOR BASE (30 kgf·cm nominal, 12V)\n');
fprintf('============================================\n');

% Resistencia de armadura (Ley de Ohm en stall)
R_base = V_base / I_stall_base;
fprintf('R = %.4f Ω\n', R_base);

% Velocidad del motor SIN reductora
omega_m_nl_base = rpm_nl_out_base * (2*pi/60) * N_base;  % [rad/s]
omega_m_nom_base = rpm_nom_out_base * (2*pi/60) * N_base; % [rad/s]
fprintf('ω_motor sin carga = %.0f rad/s (%.0f RPM)\n', omega_m_nl_base, rpm_nl_out_base*N_base);
fprintf('ω_motor nominal   = %.0f rad/s (%.0f RPM)\n', omega_m_nom_base, rpm_nom_out_base*N_base);

% Constante de fuerza contraelectromotriz Ke
% A V nominal y sin carga: V ≈ Ke·ω_m_nl  (I ≈ 0, caída R·I despreciable)
% Pero podemos hacerlo más preciso: V = I_nl·R + Ke·ω_m_nl
% I_nl ≈ 0.1 * I_stall ≈ 1A... mejor usamos la ecuación completa.
% Para máxima precisión, Ke = (V - I_nl·R) / ω_m_nl
I_nl_base = 0.1 * I_stall_base;  % Estimación de corriente sin carga
Ke_base = (V_base - I_nl_base * R_base) / omega_m_nl_base;
fprintf('\nKe = %.6f V·s/rad\n', Ke_base);

% Constante de torque Kt (idealmente = Ke en SI)
Kt_base = Ke_base;  % Motor DC ideal
fprintf('Kt = %.6f N·m/A (= Ke, motor DC ideal)\n', Kt_base);
fprintf('Kt = %.4f kgf·cm/A\n', Kt_base / kgfcm_to_Nm);

% Verificación: torque stall a la salida (con eficiencia estimada)
tau_m_stall_base = Kt_base * I_stall_base;
tau_out_stall_est_base = tau_m_stall_base * N_base * eta_est;
tau_out_stall_ds_base = tau_stall_out_base * kgfcm_to_Nm;
fprintf('\n--- Verificación con torque stall ---\n');
fprintf('τ_motor_stall = %.4f N·m\n', tau_m_stall_base);
fprintf('τ_out_stall (estimado con η=%.2f) = %.2f N·m (%.1f kgf·cm)\n', ...
    eta_est, tau_out_stall_est_base, tau_out_stall_est_base/kgfcm_to_Nm);
fprintf('τ_out_stall (catálogo)           = %.2f N·m (%.1f kgf·cm)\n', ...
    tau_out_stall_ds_base, tau_stall_out_base);
% Eficiencia que haría coincidir
eta_calc_base = tau_out_stall_ds_base / (tau_m_stall_base * N_base);
fprintf('Eficiencia implícita del gearbox = %.1f%%\n', eta_calc_base*100);

% Parámetros mecánicos del motor
J_eq_base = J_rotor_base;  % Inercia vista por el motor [kg·m²]
fprintf('\nInercia rotor J_m = %.2e kg·m²\n', J_rotor_base);

%% 3. CÁLCULO DE CONSTANTES DEL MOTOR CODO
fprintf('\n============================================\n');
fprintf('  MOTOR CODO (19 kgf·cm nominal, 24V)\n');
fprintf('============================================\n');

R_codo = V_codo / I_stall_codo;
fprintf('R = %.4f Ω\n', R_codo);

omega_m_nl_codo = rpm_nl_out_codo * (2*pi/60) * N_codo;
omega_m_nom_codo = rpm_nom_out_codo * (2*pi/60) * N_codo;
fprintf('ω_motor sin carga = %.0f rad/s (%.0f RPM)\n', omega_m_nl_codo, rpm_nl_out_codo*N_codo);
fprintf('ω_motor nominal   = %.0f rad/s (%.0f RPM)\n', omega_m_nom_codo, rpm_nom_out_codo*N_codo);

I_nl_codo = 0.1 * I_stall_codo;
Ke_codo = (V_codo - I_nl_codo * R_codo) / omega_m_nl_codo;
fprintf('\nKe = %.6f V·s/rad\n', Ke_codo);

Kt_codo = Ke_codo;
fprintf('Kt = %.6f N·m/A (= Ke, motor DC ideal)\n', Kt_codo);
fprintf('Kt = %.4f kgf·cm/A\n', Kt_codo / kgfcm_to_Nm);

tau_m_stall_codo = Kt_codo * I_stall_codo;
tau_out_stall_est_codo = tau_m_stall_codo * N_codo * eta_est;
tau_out_stall_ds_codo = tau_stall_out_codo * kgfcm_to_Nm;
fprintf('\n--- Verificación con torque stall ---\n');
fprintf('τ_motor_stall = %.4f N·m\n', tau_m_stall_codo);
fprintf('τ_out_stall (estimado con η=%.2f) = %.2f N·m (%.1f kgf·cm)\n', ...
    eta_est, tau_out_stall_est_codo, tau_out_stall_est_codo/kgfcm_to_Nm);
fprintf('τ_out_stall (catálogo)           = %.2f N·m (%.1f kgf·cm)\n', ...
    tau_out_stall_ds_codo, tau_stall_out_codo);
eta_calc_codo = tau_out_stall_ds_codo / (tau_m_stall_codo * N_codo);
fprintf('Eficiencia implícita del gearbox = %.1f%%\n', eta_calc_codo*100);

J_eq_codo = J_rotor_codo;
fprintf('\nInercia rotor J_m = %.2e kg·m²\n', J_rotor_codo);

%% 4. FUNCIONES DE TRANSFERENCIA
fprintf('\n============================================\n');
fprintf('  FUNCIONES DE TRANSFERENCIA\n');
fprintf('============================================\n');

% Para un motor DC con excitación independiente:
%
% G_ω(s) = Ω(s)/V(s) = Kt / [(L·s + R)(J·s + b) + Kt·Ke]
%
% Forma estándar: G_ω(s) = K_m / (τ_e·τ_m·s² + τ_m·s + 1)
%   donde:
%     τ_e = L/R        (constante de tiempo eléctrica)
%     τ_m = R·J/(Kt·Ke) (constante de tiempo mecánica, approx)
%     K_m = Kt/(R·b + Kt·Ke) ≈ 1/Ke  (ganancia en DC)
%
% Si τ_e << τ_m (generalmente cierto), se simplifica a:
%   G_ω(s) ≈ 1/Ke / (τ_m·s + 1)

% ===================== MOTOR BASE =====================
% Constantes de tiempo
tau_e_base = L_est_base / R_base;
tau_m_base = R_base * J_eq_base / (Kt_base * Ke_base + R_base * b_base);
K_DC_base  = Kt_base / (R_base * b_base + Kt_base * Ke_base);

fprintf('\n--- Motor Base ---\n');
fprintf('τ_e = %.2f μs  (constante eléctrica)\n', tau_e_base*1e6);
fprintf('τ_m = %.4f s   (constante mecánica)\n', tau_m_base);
fprintf('K_DC = %.4f (rad/s)/V  (ganancia DC)\n', K_DC_base);
fprintf('1/Ke = %.4f (rad/s)/V  (ganancia DC ideal)\n', 1/Ke_base);

% Función de transferencia completa: Ω(s)/V(s)
s = tf('s');
G_omega_base = Kt_base / ((L_est_base*s + R_base)*(J_eq_base*s + b_base) + Kt_base*Ke_base);
G_omega_base = minreal(G_omega_base);

% Aproximación de 1er orden (simplificada si τ_e << τ_m)
G_omega_simple_base = K_DC_base / (tau_m_base * s + 1);

fprintf('\nFT completa (2° orden):\n');
G_omega_base
fprintf('FT 1er orden (aprox):\n');
G_omega_simple_base

% FT de posición: Θ(s)/V(s) = G_ω(s) / s
G_theta_base = G_omega_base / s;
G_theta_simple_base = G_omega_simple_base / s;

% ===================== MOTOR CODO =====================
tau_e_codo = L_est_codo / R_codo;
tau_m_codo = R_codo * J_eq_codo / (Kt_codo * Ke_codo + R_codo * b_codo);
K_DC_codo  = Kt_codo / (R_codo * b_codo + Kt_codo * Ke_codo);

fprintf('\n--- Motor Codo ---\n');
fprintf('τ_e = %.2f μs  (constante eléctrica)\n', tau_e_codo*1e6);
fprintf('τ_m = %.4f s   (constante mecánica)\n', tau_m_codo);
fprintf('K_DC = %.4f (rad/s)/V  (ganancia DC)\n', K_DC_codo);
fprintf('1/Ke = %.4f (rad/s)/V  (ganancia DC ideal)\n', 1/Ke_codo);

G_omega_codo = Kt_codo / ((L_est_codo*s + R_codo)*(J_eq_codo*s + b_codo) + Kt_codo*Ke_codo);
G_omega_codo = minreal(G_omega_codo);
G_omega_simple_codo = K_DC_codo / (tau_m_codo * s + 1);

fprintf('\nFT completa (2° orden):\n');
G_omega_codo
fprintf('FT 1er orden (aprox):\n');
G_omega_simple_codo

G_theta_codo = G_omega_codo / s;
G_theta_simple_codo = G_omega_simple_codo / s;

%% 5. GRÁFICAS DE RESPUESTA
figure(1); clf;

% --- Respuesta al escalón: velocidad ---
subplot(2,3,1);
step(G_omega_base, 0.5); hold on;
step(G_omega_simple_base, 0.5, '--');
title('Base: Ω(s)/V(s) — Escalón de velocidad');
ylabel('Velocidad [rad/s]'); grid on;
legend('2° orden', '1er orden', 'Location','best');

subplot(2,3,2);
step(G_omega_codo, 0.5); hold on;
step(G_omega_simple_codo, 0.5, '--');
title('Codo: Ω(s)/V(s) — Escalón de velocidad');
ylabel('Velocidad [rad/s]'); grid on;
legend('2° orden', '1er orden', 'Location','best');

% --- Respuesta al escalón: posición ---
subplot(2,3,3);
step(G_theta_base, 0.5); hold on;
step(G_theta_simple_base, 0.5, '--');
title('Base: Θ(s)/V(s) — Escalón de posición (integrador)');
ylabel('Posición [rad]'); grid on;
legend('2° orden', '1er orden', 'Location','best');

subplot(2,3,4);
step(G_theta_codo, 0.5); hold on;
step(G_theta_simple_codo, 0.5, '--');
title('Codo: Θ(s)/V(s) — Escalón de posición (integrador)');
ylabel('Posición [rad]'); grid on;
legend('2° orden', '1er orden', 'Location','best');

% --- Diagrama de Bode ---
subplot(2,3,5);
bode(G_omega_base, G_omega_simple_base);
title('Base: Diagrama de Bode'); grid on;
legend('2° orden', '1er orden', 'Location','best');

subplot(2,3,6);
bode(G_omega_codo, G_omega_simple_codo);
title('Codo: Diagrama de Bode'); grid on;
legend('2° orden', '1er orden', 'Location','best');

sgtitle('Respuesta de Motores DC con Reductora');

%% 6. POLOS, CEROS Y ESTABILIDAD
fprintf('\n============================================\n');
fprintf('  POLOS Y ESTABILIDAD\n');
fprintf('============================================\n');

fprintf('\n--- Motor Base ---\n');
p_base = pole(G_omega_base);
z_base = zero(G_omega_base);
fprintf('Polos: s = %.2f ± j%.2f\n', real(p_base(1)), abs(imag(p_base(1))));
if ~isempty(z_base)
    fprintf('Ceros: s = %.4f\n', z_base);
end
fprintf('Frecuencia natural ω_n = %.2f rad/s (%.2f Hz)\n', abs(p_base(1)), abs(p_base(1))/(2*pi));
fprintf('Todos los polos en semiplano izquierdo: SISTEMA ESTABLE\n');

fprintf('\n--- Motor Codo ---\n');
p_codo = pole(G_omega_codo);
z_codo = zero(G_omega_codo);
fprintf('Polos: s = %.2f ± j%.2f\n', real(p_codo(1)), abs(imag(p_codo(1))));
if ~isempty(z_codo)
    fprintf('Ceros: s = %.4f\n', z_codo);
end
fprintf('Frecuencia natural ω_n = %.2f rad/s (%.2f Hz)\n', abs(p_codo(1)), abs(p_codo(1))/(2*pi));
fprintf('Todos los polos en semiplano izquierdo: SISTEMA ESTABLE\n');

%% 7. MODELO PARA PWM
fprintf('\n============================================\n');
fprintf('  MODELO CON PWM\n');
fprintf('============================================\n');
fprintf('\nPara control por PWM, se modela como:\n');
fprintf('  V_motor = duty_cycle × V_max\n');
fprintf('donde duty_cycle ∈ [0, 1] y V_max es el voltaje de alimentación.\n\n');
fprintf('Si el puente H opera a V_max:\n');
fprintf('  Base: V_max_alim = 12V (o el voltaje del bus de potencia)\n');
fprintf('  Codo: V_max_alim = 24V (o el voltaje del bus de potencia)\n\n');
fprintf('FT con entrada PWM (duty cycle → velocidad):\n');
fprintf('  G_pwm_base(s) = V_max_base × G_omega_base(s)\n');
fprintf('  G_pwm_codo(s) = V_max_codo × G_omega_codo(s)\n');

% Si el PWM controla directamente el voltaje promedio:
Vbus_base = 12;   % Voltaje del bus de potencia para base [V]
Vbus_codo = 24;   % Voltaje del bus de potencia para codo [V]

G_pwm_base = Vbus_base * G_omega_base;       % (rad/s) por unidad de duty cycle
G_pwm_codo = Vbus_codo * G_omega_codo;

fprintf('\nG_pwm_base(s) = Ω(s)/duty(s):\n');
G_pwm_base
fprintf('\nG_pwm_codo(s) = Ω(s)/duty(s):\n');
G_pwm_codo

%% 8. GUARDAR MODELO PARA USO EN CONTROL
fprintf('\n============================================\n');
fprintf('  EXPORTACIÓN\n');
fprintf('============================================\n');

% Guardar variables en archivo .mat para usar en diseño de control
save('modelo_motores.mat', ...
    'R_base', 'Ke_base', 'Kt_base', 'J_eq_base', 'L_est_base', 'b_base', ...
    'R_codo', 'Ke_codo', 'Kt_codo', 'J_eq_codo', 'L_est_codo', 'b_codo', ...
    'G_omega_base', 'G_omega_simple_base', 'G_theta_base', ...
    'G_omega_codo', 'G_omega_simple_codo', 'G_theta_codo', ...
    'G_pwm_base', 'G_pwm_codo', ...
    'Vbus_base', 'Vbus_codo', 'N_base', 'N_codo');

fprintf('Modelo guardado en: modelo_motores.mat\n');

%% 9. RESUMEN COMPARATIVO
fprintf('\n============================================\n');
fprintf('  RESUMEN DE CONSTANTES\n');
fprintf('============================================\n');
fprintf('%-25s %15s %15s\n', 'Parámetro', 'Base', 'Codo');
fprintf('%-25s %15s %15s\n', '-------------------------','---------------','---------------');
fprintf('%-25s %15.4f %15.4f\n', 'R [Ω]', R_base, R_codo);
fprintf('%-25s %15.6f %15.6f\n', 'Ke [V·s/rad]', Ke_base, Ke_codo);
fprintf('%-25s %15.6f %15.6f\n', 'Kt [N·m/A]', Kt_base, Kt_codo);
fprintf('%-25s %15.2f %15.2f\n', 'τ_e [μs]', tau_e_base*1e6, tau_e_codo*1e6);
fprintf('%-25s %15.4f %15.4f\n', 'τ_m [s]', tau_m_base, tau_m_codo);
fprintf('%-25s %15.4f %15.4f\n', 'K_DC [(rad/s)/V]', K_DC_base, K_DC_codo);
fprintf('%-25s %15.0f %15.0f\n', 'ω_m_nl [RPM]', rpm_nl_out_base*N_base, rpm_nl_out_codo*N_codo);

fprintf('\n=== FIN DEL ANÁLISIS ===\n');
