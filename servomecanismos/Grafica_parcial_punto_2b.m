% Script MATLAB - Respuesta a Escalon para Sistema Servo (Punto 2b)
% Servomecanismos v2025-2s
% Basado en diagrama de bloques: PID + FF + Motor DC + STPM + Feedback

clc; clear; close all;

%% 1. Definición de Parámetros del Sistema
% Motor DC (Valores típicos para un servomotor pequeño)
R_a = 2.0;      % Resistencia de armadura [Ohm]
K_t = 0.05;     % Constante de torque [Nm/A]
K_e = 0.05;     % Constante de FEM [V/(rad/s)]
J_m = 1e-4;     % Inercia del motor [kg*m^2]
b_m = 1e-5;     % Coeficiente de fricción viscosa [Nm/(rad/s)]

% STPM (Husillo de bolas)
p = 0.01;       % Paso del husillo [m/rev]
K_husillo = p / (2*pi); % Ganancia de conversión rotación->lineal [m/rad]

% Controlador PID (Sintonizado para respuesta estable)
Kp = 150;       % Ganancia Proporcional
Ki = 50;        % Ganancia Integral
Kd = 0.5;       % Ganancia Derivativa

% Controlador Feedforward (Simplificado para este ejemplo)
FF_gain = 0.5;  % Porcentaje de acción directa

% Driver y Sensor (Ganancias idealizadas)
K_driver = 1;
H_sensor = 1;   % Retroalimentación unitaria

%% 2. Construcción de la Función de Transferencia

% --- Planta del Motor (Velocidad angular / Voltaje) ---
% Ecuación: (J_m*s + b_m)*omega = (K_t * I_a)
% I_a = (V_in - K_e*omega) / R_a
% Despejando: G_motor(s) = omega(s)/V_in(s) = K_t / (R_a*J_m*s + R_a*b_m + K_t*K_e)
s = tf('s');
G_motor = K_t / (R_a * J_m * s + R_a * b_m + K_t * K_e);

% --- Planta Mecánica (Posición angular / Velocidad angular) ---
G_integrador = 1 / s;

% --- Planta Husillo (Posición lineal / Posición angular) ---
% G_husillo = p / (2*pi)
G_husillo = K_husillo;

% --- Planta Total (Sin Controlador) ---
% Desde Voltaje de entrada del motor hasta Posición lineal de salida
G_planta = G_motor * G_integrador * G_husillo;

% --- Controlador Total (PID + Feedforward) ---
C_pid = Kp + Ki/s + Kd*s;
C_total = C_pid + FF_gain; % Superposición de FF y PID como en el diagrama

% --- Sistema en Lazo Cerrado ---
% T(s) = (C_total * K_driver * G_planta) / (1 + (C_total * K_driver * G_planta * H_sensor))
G_lazo_abierto = C_total * K_driver * G_planta;
T = feedback(G_lazo_abierto, H_sensor);

%% 3. Simulación y Gráfica
% Entrada Escalón Unitario de Posición (Ref = 1 metro)
t = 0:0.001:2; % Tiempo de simulación 2 segundos
[y, t] = step(T, t);

% Crear figura
figure('Name', 'Respuesta a Escalón - Sistema Servo', 'Color', 'w');
plot(t, y, 'LineWidth', 2);
grid on;
hold on;

%% 4. Estimación de Parámetros (Anotaciones en la gráfica)

% -- Estimación de Ganancia (Valor final en régimen permanente) --
y_final = y(end); % Valor asintótico
yline(y_final, '--r', 'LineWidth', 1.5);
text(1.2, y_final*0.95, sprintf('Ganancia K \\approx %.2f', y_final), 'FontSize', 10, 'Color', 'r');

% -- Estimación de Constante de tiempo (63.2% del valor final) --
y_63 = 0.632 * y_final;
% Encontrar el índice donde se cruza el 63.2%
idx = find(y >= y_63, 1, 'first');
tau = t(idx); % Tiempo estimado

yline(y_63, '--g', 'LineWidth', 1.5);
xline(tau, '--g', 'LineWidth', 1.5);
plot(tau, y_63, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
text(tau + 0.05, y_63 - 0.05, sprintf('\\tau \\approx %.3f s', tau), 'FontSize', 10, 'Color', 'g');

% -- Anotaciones adicionales --
title('Respuesta a Escalón Unitario de Posición');
xlabel('Tiempo (s)');
ylabel('Posición (m)');
legend('Respuesta del sistema', 'Valor Final (Ganancia)', '63.2% (Constante de tiempo)', 'Location', 'southeast');
ylim([0, 1.2]);

% Mostrar resultados en consola
fprintf('--- RESULTADOS DE SIMULACIÓN ---\n');
fprintf('Ganancia del sistema (K): %.4f\n', y_final);
fprintf('Constante de tiempo (tau): %.4f s\n', tau);
fprintf('Tiempo de establecimiento (2%%): %.4f s\n', find(y >= 0.98*y_final, 1, 'first')*0.001);

disp('Gráfica generada exitosamente.');