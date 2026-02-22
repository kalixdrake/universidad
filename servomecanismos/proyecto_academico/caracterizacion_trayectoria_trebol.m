%% CARACTERIZACIÓN, CINEMÁTICA Y DINÁMICA DE ROBOT 2R (OPTIMIZADO PARA ESP32)
clear; clc; close all;

%% 1. Parámetros del Proyecto 
L = 0.20;           % Lado mínimo de 20cm
n = 7;              % Número de hojas
v_const = 0.1;     % Velocidad constante deseada (1 cm/s)
S = 1.25;           % Escala inicial
beta = deg2rad(30);  % Rotación inicial

% Longitud de eslabones y masas
l1 = 0.18; m_link1 = 0.700;
l2 = 0.18; m_link2 = 0.600;
m_motor2 = 0.300; m_tip = 0.050;
g = 9.81;

%% 2. Generación de la Trayectoria Base del Trébol
r0 = (L/2) * 0.75;  
A = (L/2) * 0.25;   

phi_raw = linspace(0, 2*pi, 1000); 
r_raw = r0 + A*cos(n * phi_raw);
x_raw = r_raw .* cos(phi_raw);
y_raw = r_raw .* sin(phi_raw);

% Parametrización para velocidad constante
dx = diff(x_raw); dy = diff(y_raw);
ds = sqrt(dx.^2 + dy.^2);
s_acumulada = [0, cumsum(ds)];
T_total = s_acumulada(end) / v_const;
dt = 0.01; % Paso de tiempo de simulación
t = 0:dt:T_total;

s_t = v_const * t; 
phi_t = interp1(s_acumulada, phi_raw, s_t, 'pchip', 'extrap');

r_t = r0 + A*cos(n * phi_t);
x_traj = S * r_t .* cos(phi_t + beta);
y_traj = S * r_t .* sin(phi_t + beta);
x_traj(end) = x_traj(1); y_traj(end) = y_traj(1);

% Desplazamiento del Origen al espacio de trabajo seguro
X_centro = 0.1;  Y_centro = 0.19; 
x_traj_desp = x_traj + X_centro;
y_traj_desp = y_traj + Y_centro;

%% 3. Optimización de la Posición "Home" (Codo Arriba)
disp('1. Calculando posición Home óptima...');
th1_home = pi/2; % Eslabón 1 vertical
Y_codo = l1 * sin(th1_home);

th2_test = linspace(-pi, 0, 180); % Solo codo arriba
Y_efector_test = Y_codo + l2 * sin(th1_home + th2_test);

posiciones_validas = (Y_efector_test >= 0.02) & (Y_efector_test <= Y_centro);

% Calcular torques estáticos aproximados para elegir el menor esfuerzo
m2_t = m_link2 + m_tip;
lc2_t = (m_link2 * (l2/2) + m_tip * l2) / m2_t;
G1_test = m2_t * lc2_t * g * cos(th1_home + th2_test);
G2_test = m2_t * lc2_t * g * cos(th1_home + th2_test);
Torque_Estacionario = abs(G1_test) + abs(G2_test);
Torque_Estacionario(~posiciones_validas) = inf;

[~, idx] = min(Torque_Estacionario);
th2_home = th2_test(idx);

x_inicio = l1 * cos(th1_home) + l2 * cos(th1_home + th2_home);
y_inicio = l1 * sin(th1_home) + l2 * sin(th1_home + th2_home);

%% 4. Reordenar Trébol para Iniciar en el Punto Más Cercano
disp('2. Alineando trayectoria con el Home...');
distancias_al_home = sqrt((x_traj_desp - x_inicio).^2 + (y_traj_desp - y_inicio).^2);
[~, idx_closest] = min(distancias_al_home);

x_temp = x_traj_desp(1:end-1); y_temp = y_traj_desp(1:end-1);
x_traj_opt = [x_temp(idx_closest:end), x_temp(1:idx_closest-1)];
y_traj_opt = [y_temp(idx_closest:end), y_temp(1:idx_closest-1)];

x_traj_desp = [x_traj_opt, x_traj_opt(1)];
y_traj_desp = [y_traj_opt, y_traj_opt(1)];

%% 5. Fase de Aproximación: Línea Recta con Arranque Suave y Pausa
disp('3. Generando aproximación inicial suave...');
T_pausa = 1.0; 
t_pausa = 0:dt:T_pausa;
x_pausa = x_inicio * ones(1, length(t_pausa));
y_pausa = y_inicio * ones(1, length(t_pausa));
x_pausa(end) = []; y_pausa(end) = [];

dist_aprox = sqrt((x_traj_desp(1) - x_inicio)^2 + (y_traj_desp(1) - y_inicio)^2);
T_aprox = dist_aprox / v_const; 
t_aprox = 0:dt:T_aprox;

tau = t_aprox / T_aprox; 
s_tau = 3*tau.^2 - 2*tau.^3; % Perfil de velocidad en "S"

x_aprox = x_inicio + (x_traj_desp(1) - x_inicio) * s_tau;
y_aprox = y_inicio + (y_traj_desp(1) - y_inicio) * s_tau;
x_aprox(end) = []; y_aprox(end) = [];

x_total = [x_pausa, x_aprox, x_traj_desp];
y_total = [y_pausa, y_aprox, y_traj_desp];
num_puntos = length(x_total);

%% 6. Cinemática Inversa Vectorizada (Codo Arriba Garantizado)
disp('4. Calculando Cinemática Inversa...');
theta1_raw = zeros(1, num_puntos);
theta2_raw = zeros(1, num_puntos);

for i = 1:num_puntos
    x = x_total(i); y = y_total(i);
    D = (x^2 + y^2 - l1^2 - l2^2) / (2 * l1 * l2);
    if D > 1, D = 1; elseif D < -1, D = -1; end
    
    th2 = atan2(-sqrt(1 - D^2), D);  
    th1 = atan2(y, x) - atan2(l2 * sin(th2), l1 + l2 * cos(th2));
    
    theta1_raw(i) = th1; theta2_raw(i) = th2;
end
theta1_raw = unwrap(theta1_raw);
theta2_raw = unwrap(theta2_raw);

%% 7. Filtrado Savitzky-Golay (Aplanado de Aceleraciones y Jerk)
disp('5. Optimizando perfiles para motores (Filtrado Savitzky-Golay)...');
orden_pol = 3;       
window_size = 5;    % Ajusta este valor (impar) si quieres más o menos suavizado

theta1 = sgolayfilt(theta1_raw, orden_pol, window_size);
theta2 = sgolayfilt(theta2_raw, orden_pol, window_size);

dtheta1 = sgolayfilt(gradient(theta1, dt), orden_pol, window_size);
dtheta2 = sgolayfilt(gradient(theta2, dt), orden_pol, window_size);

ddtheta1 = sgolayfilt(gradient(dtheta1, dt), orden_pol, window_size);
ddtheta2 = sgolayfilt(gradient(dtheta2, dt), orden_pol, window_size);

t_total_vec = (0:num_puntos-1) * dt;

%% 8. Visualización de la Deformación del Filtrado
x_ideal = l1 * cos(theta1_raw) + l2 * cos(theta1_raw + theta2_raw);
y_ideal = l1 * sin(theta1_raw) + l2 * sin(theta1_raw + theta2_raw);

x_real = l1 * cos(theta1) + l2 * cos(theta1 + theta2);
y_real = l1 * sin(theta1) + l2 * sin(theta1 + theta2);

figure(1); clf; hold on; axis equal; grid on;
plot(x_ideal, y_ideal, 'Color', [0.8 0.8 0.8], 'LineWidth', 3, 'DisplayName', 'Trébol Matemático Ideal');
plot(x_real, y_real, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Trébol Suavizado (Motores)');
title('Comparación: Trayectoria Ideal vs Suavizada');
legend('Location', 'best');

%% 9. Dinámica Inversa y Cálculo de Torques
disp('6. Calculando Dinámica Inversa...');
m1_t = m_link1 + m_motor2;
lc1_t = (m_link1 * (l1/2) + m_motor2 * l1) / m1_t; 
I1_t = (1/12)*m_link1*l1^2 + m_link1*(lc1_t - l1/2)^2 + m_motor2*(l1 - lc1_t)^2; 

I2_t = (1/12)*m_link2*l2^2 + m_link2*(lc2_t - l2/2)^2 + m_tip*(l2 - lc2_t)^2;

tau1 = zeros(1, num_puntos); tau2 = zeros(1, num_puntos);

for i = 1:num_puntos
    th1 = theta1(i); th2 = theta2(i);
    d1 = dtheta1(i); d2 = dtheta2(i);
    dd1 = ddtheta1(i); dd2 = ddtheta2(i);
    
    M11 = m1_t*lc1_t^2 + m2_t*(l1^2 + lc2_t^2 + 2*l1*lc2_t*cos(th2)) + I1_t + I2_t;
    M12 = m2_t*(lc2_t^2 + l1*lc2_t*cos(th2)) + I2_t;
    M22 = m2_t*lc2_t^2 + I2_t;
    
    h_coriolis = -m2_t*l1*lc2_t*sin(th2);
    C1 = h_coriolis*(2*d1*d2 + d2^2);
    C2 = -h_coriolis*d1^2;
    
    G1_dyn = (m1_t*lc1_t + m2_t*l1)*g*cos(th1) + m2_t*lc2_t*g*cos(th1 + th2);
    G2_dyn = m2_t*lc2_t*g*cos(th1 + th2);
    
    tau1(i) = M11*dd1 + M12*dd2 + C1 + G1_dyn;
    tau2(i) = M12*dd1 + M22*dd2 + C2 + G2_dyn;
end

tau1_kgfcm = tau1 * (100 / 9.81);
tau2_kgfcm = tau2 * (100 / 9.81);


%% 10. Gráficas Cinemáticas y de Torque
figure(2); clf;
subplot(3,1,1); plot(t_total_vec, dtheta1, 'b', t_total_vec, dtheta2, 'r'); title('Velocidades (rad/s)'); grid on;
subplot(3,1,2); plot(t_total_vec, ddtheta1, 'b', t_total_vec, ddtheta2, 'r'); title('Aceleraciones Suavizadas (rad/s^2)'); grid on;
subplot(3,1,3); plot(t_total_vec, tau1_kgfcm, 'b', t_total_vec, tau2_kgfcm, 'r'); title('Torques (kgf·cm)'); xlabel('Tiempo [s]'); grid on;

%% 11. Animación del Robot
figure(3); clf; hold on; axis equal; grid on; axis([-0.25 0.6 -0.05 0.6]);
plot([-0.25, 0.6], [0, 0], 'k-', 'LineWidth', 3); % Mesa
plot(x_total, y_total, 'color', [0.8 0.8 0.8], 'LineWidth', 1.5); % Ruta

h_robot = plot([0, l1*cos(theta1(1)), x_real(1)], [0, l1*sin(theta1(1)), y_real(1)], '-o', 'LineWidth', 4, 'Color', '#D95319', 'MarkerFaceColor', '#0072BD');
h_rastro = plot(x_real(1), y_real(1), 'b-', 'LineWidth', 2);

x_rastro = zeros(1, num_puntos); y_rastro = zeros(1, num_puntos);
paso = 15;

disp('7. Reproduciendo animación...');
for i = 1:paso:num_puntos
    x_c = l1*cos(theta1(i)); y_c = l1*sin(theta1(i));
    x_rastro(1:i) = x_real(1:i); y_rastro(1:i) = y_real(1:i);
    
    set(h_robot, 'XData', [0, x_c, x_real(i)], 'YData', [0, y_c, y_real(i)]);
    set(h_rastro, 'XData', x_rastro(1:i), 'YData', y_rastro(1:i));
    drawnow;
end
set(h_robot, 'XData', [0, l1*cos(theta1(end)), x_real(end)], 'YData', [0, l1*sin(theta1(end)), y_real(end)]);
set(h_rastro, 'XData', x_real, 'YData', y_real); drawnow;
disp('Proceso finalizado.');

%% 12. BARRIDO DE PARÁMETROS: BÚSQUEDA DEL CASO CRÍTICO (WORST-CASE SCENARIO)
disp('Iniciando análisis del Caso Crítico. Esto puede tardar unos segundos...');

% 1. Definir los vectores de prueba (resolución ajustable)
vec_S = linspace(0.75, 1.25, 10);        % 10 escalas
vec_n = 4:7;                            % 4 tipos de trébol
vec_v = linspace(0.01, 0.1, 10);         % 10 velocidades
vec_beta = deg2rad(linspace(-45, 45, 90)); % 90 rotaciones

% 2. Variables para rastrear los máximos históricos
max_global_tau1 = 0; max_global_tau2 = 0;
max_global_acc1 = 0; max_global_acc2 = 0;
params_criticos = struct('S', 0, 'n', 0, 'v', 0, 'beta', 0);

total_iter = length(vec_S) * length(vec_n) * length(vec_v) * length(vec_beta);
iter_actual = 0;

for S_i = vec_S
    for n_i = vec_n
        for v_i = vec_v
            for beta_i = vec_beta
                iter_actual = iter_actual + 1;
                
                % --- A. Generación de Trayectoria Rápida ---
                phi_raw = linspace(0, 2*pi, 500); % Menos puntos para iterar rápido
                r_raw = r0 + A*cos(n_i * phi_raw);
                x_raw = r_raw .* cos(phi_raw);
                y_raw = r_raw .* sin(phi_raw);
                
                dx = diff(x_raw); dy = diff(y_raw);
                ds = sqrt(dx.^2 + dy.^2);
                s_acum = [0, cumsum(ds)];
                
                T_tot = s_acum(end) / v_i;
                dt_sim = 0.01; 
                t_sim = 0:dt_sim:T_tot;
                s_t = v_i * t_sim; 
                
                phi_t = interp1(s_acum, phi_raw, s_t, 'pchip', 'extrap');
                r_t = r0 + A*cos(n_i * phi_t);
                
                x_traj = S_i * r_t .* cos(phi_t + beta_i) + X_centro;
                y_traj = S_i * r_t .* sin(phi_t + beta_i) + Y_centro;
                
                % --- B. Comprobación de Espacio de Trabajo ---
                D_test = (x_traj.^2 + y_traj.^2 - l1^2 - l2^2) / (2 * l1 * l2);
                if any(D_test > 1) || any(D_test < -1)
                    continue; % Si la trayectoria sale del alcance, la saltamos
                end
                
                % --- C. Cinemática Inversa ---
                th2_iter = atan2(-sqrt(1 - D_test.^2), D_test);  
                th1_iter = atan2(y_traj, x_traj) - atan2(l2 * sin(th2_iter), l1 + l2 * cos(th2_iter));
                
                th1_iter = unwrap(th1_iter); th2_iter = unwrap(th2_iter);
                
                % --- D. Filtrado de Velocidad y Aceleración ---
                th1_filt = sgolayfilt(th1_iter, 3, 21);
                th2_filt = sgolayfilt(th2_iter, 3, 21);
                
                d1_iter = sgolayfilt(gradient(th1_filt, dt_sim), 3, 21);
                d2_iter = sgolayfilt(gradient(th2_filt, dt_sim), 3, 21);
                
                dd1_iter = sgolayfilt(gradient(d1_iter, dt_sim), 3, 21);
                dd2_iter = sgolayfilt(gradient(d2_iter, dt_sim), 3, 21);
                
                % --- E. Dinámica Inversa Rápida ---
                M11 = m1_t*lc1_t^2 + m2_t*(l1^2 + lc2_t^2 + 2*l1*lc2_t*cos(th2_filt)) + I1_t + I2_t;
                M12 = m2_t*(lc2_t^2 + l1*lc2_t*cos(th2_filt)) + I2_t;
                M22 = m2_t*lc2_t^2 + I2_t;
                
                h_cor = -m2_t*l1*lc2_t*sin(th2_filt);
                C1 = h_cor.*(2.*d1_iter.*d2_iter + d2_iter.^2);
                C2 = -h_cor.*d1_iter.^2;
                
                G1_dyn = (m1_t*lc1_t + m2_t*l1)*g*cos(th1_filt) + m2_t*lc2_t*g*cos(th1_filt + th2_filt);
                G2_dyn = m2_t*lc2_t*g*cos(th1_filt + th2_filt);
                
                tau1_iter = M11.*dd1_iter + M12.*dd2_iter + C1 + G1_dyn;
                tau2_iter = M12.*dd1_iter + M22.*dd2_iter + C2 + G2_dyn;
                
                tau1_kgfcm_iter = tau1_iter * (100 / 9.81);
                tau2_kgfcm_iter = tau2_iter * (100 / 9.81);
                
                % --- F. Evaluar y Actualizar Máximos ---
                max_t1_actual = max(abs(tau1_kgfcm_iter));
                max_t2_actual = max(abs(tau2_kgfcm_iter));
                
                if max_t1_actual > max_global_tau1 || max_t2_actual > max_global_tau2
                    max_global_tau1 = max(max_global_tau1, max_t1_actual);
                    max_global_tau2 = max(max_global_tau2, max_t2_actual);
                    max_global_acc1 = max(max_global_acc1, max(abs(dd1_iter)));
                    max_global_acc2 = max(max_global_acc2, max(abs(dd2_iter)));
                    
                    params_criticos.S = S_i;
                    params_criticos.n = n_i;
                    params_criticos.v = v_i;
                    params_criticos.beta = rad2deg(beta_i);
                end
            end
        end
    end
end

%% 13. REPORTE DEL CASO CRÍTICO
fprintf('\n=======================================================\n');
fprintf('   REPORTE DE DIMENSIONAMIENTO: CASO CRÍTICO ABSOLUTO\n');
fprintf('=======================================================\n');
fprintf('Configuración más exigente encontrada:\n');
fprintf('  - Escala (S): %.2f\n', params_criticos.S);
fprintf('  - Número de Hojas (n): %d\n', params_criticos.n);
fprintf('  - Velocidad (v): %.3f m/s\n', params_criticos.v);
fprintf('  - Rotación Inicial: %.1f°\n\n', params_criticos.beta);

fprintf('ESPECIFICACIONES MÍNIMAS PARA COMPRA DE MOTORES:\n');
fprintf('  MOTOR 1 (BASE):\n');
fprintf('    Aceleración Pico: %.2f rad/s^2\n', max_global_acc1);
fprintf('    Torque Pico:      %.2f kgf·cm\n\n', max_global_tau1);

fprintf('  MOTOR 2 (CODO):\n');
fprintf('    Aceleración Pico: %.2f rad/s^2\n', max_global_acc2);
fprintf('    Torque Pico:      %.2f kgf·cm\n', max_global_tau2);
fprintf('=======================================================\n\n');
disp('Nota: Es recomendable añadir un factor de seguridad del 20% al 30% al torque pico al elegir el motor.');