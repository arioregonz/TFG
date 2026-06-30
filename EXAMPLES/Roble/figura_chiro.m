%
%  FIGURA DEL CHIRP PARA LA MEMORIA (Sección 3.5.1)
%  Genera una figura limpia de la señal chirp: forma temporal + espectro
%  Exportar como PNG/PDF para insertar en el TFG
%
%  Señal chirp R3α (Tabla 5 tesis): 22→50 kHz, Ts=45µs, ventana gaussiana

clear; clc; close all;

% Parámetros de la señal
DT   = 2.0e-7;   % s (paso temporal)
Ts   = 45e-6;    % s (duración)
f0   = 22e3;     % Hz (frecuencia inicial)
f1   = 50e3;     % Hz (frecuencia final)
Fc   = 36e3;     % Hz (frecuencia central)

% Vector temporal y señal
Nc = round(Ts/DT);
t  = (0:Nc-1)*DT;
mu = Ts/2; sigma = Ts/6;
f_t   = (f1-f0)*t/Ts + f0;            % frecuencia instantánea
phase = 2*pi*f_t.*t;                  % fase
win   = exp(-((t-mu).^2)/(2*sigma^2)); % ventana gaussiana
chirp = cos(phase).*win;              % señal chirp enventanada

% Espectro
NFFT = 2^nextpow2(8*Nc);
S    = abs(fft(chirp, NFFT));
S    = S(1:NFFT/2)/max(S(1:NFFT/2));   % normalizado
f    = (0:NFFT/2-1)/(NFFT*DT);

% ===== FIGURA =====
figure('Position',[100 100 900 350],'Color','w');

% Panel (a): forma temporal
subplot(1,2,1);
plot(t*1e6, chirp, 'b-', 'LineWidth', 1.2); hold on;
plot(t*1e6, win, 'r--', 'LineWidth', 1.2);
plot(t*1e6, -win, 'r--', 'LineWidth', 1.2);
grid on; box on;
xlabel('Tiempo [\mus]','FontSize',11);
ylabel('Amplitud normalizada','FontSize',11);
title('(a) Forma temporal','FontSize',11);
legend('Chirp enventanado','Ventana gaussiana','Location','southwest','FontSize',8);
xlim([0 Ts*1e6]); ylim([-1.15 1.15]);

% Panel (b): espectro
subplot(1,2,2);
plot(f/1e3, S, 'b-', 'LineWidth', 1.5); hold on;
% Marcar banda -3dB
umbral = 1/sqrt(2);
yline(umbral, 'k:', 'LineWidth', 1);
idx3 = find(S > umbral);
f_lo = f(idx3(1))/1e3; f_hi = f(idx3(end))/1e3;
xline(f0/1e3, 'g--', 'LineWidth', 1);
xline(f1/1e3, 'g--', 'LineWidth', 1);
area(f(idx3)/1e3, S(idx3), 'FaceColor', [0.7 0.85 1], 'FaceAlpha', 0.5, 'EdgeColor','none');
plot(f/1e3, S, 'b-', 'LineWidth', 1.5);
grid on; box on;
xlabel('Frecuencia [kHz]','FontSize',11);
ylabel('Magnitud normalizada','FontSize',11);
title('(b) Espectro','FontSize',11);
text(Fc/1e3, 0.5, sprintf('BW_{-3dB}\n[%.0f-%.0f] kHz', f_lo, f_hi), ...
    'HorizontalAlignment','center','FontSize',8,'BackgroundColor','w');
xlim([0 80]); ylim([0 1.1]);

% Exportar (descomentar para guardar)
% print(gcf, 'chirp_memoria.png', '-dpng', '-r300');
% print(gcf, 'chirp_memoria.pdf', '-dpdf', '-bestfit');

fprintf('Figura del chirp generada.\n');
fprintf('Banda -3dB: [%.1f - %.1f] kHz\n', f_lo, f_hi);
fprintf('Para exportar, descomenta las líneas print() al final.\n');