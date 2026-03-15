% Especificaciones
fs = 75000;
n = 4;
Rp = 1;

% Tonos del estándar P25 (Hz)
f1 = 10700; f2 = 11900; f3 = 13100; f4 = 14300;

% Media Geométrica
fcA = sqrt(f1 * f2);
fcB = sqrt(f2 * f3);
fcC = sqrt(f3 * f4);

fprintf('Frecuencias de corte:\n');
fprintf('fcA: %.2f Hz, fcB: %.2f Hz, fcC: %.2f Hz\n\n', fcA, fcB, fcC);

% 3. Transformación Bilineal

Wn_bil = [fcA fcB] / (fs/2); % Frecuencias normalizadas (0 a 1)
[b_bil, a_bil] = cheby1(n, Rp, Wn_bil, 'bandpass'); 

% Visualización
[H_bil, f] = freqz(b_bil, a_bil, 1024, fs);

% Invarianza del Impulso
% Analógico
W_analog = [2*pi*fcA, 2*pi*fcB]; % Frecuencias en rad/s
[b_s, a_s] = cheby1(n/2, Rp, W_analog, 'bandpass', 's'); % n/2 porque BPF duplica el orden

% Conversion a digital con invarianza 
[b_inv, a_inv] = impinvar(b_s, a_s, fs);

% Visualización
[H_inv, ~] = freqz(b_inv, a_inv, 1024, fs);

% 5. Plots
figure;
plot(f, 20*log10(abs(H_bil)), 'LineWidth', 1.5); hold on;
plot(f, 20*log10(abs(H_inv)), '--', 'LineWidth', 1.5);
grid on;
title('Filtro Pasabanda (11.9 kHz) - IIR Orden 4');
xlabel('Frecuencia (Hz)'); ylabel('Magnitud (dB)');
legend('Transformación Bilineal', 'Invarianza del Impulso');
axis([8000 16000 -60 5]); % Zoom en la zona de interés