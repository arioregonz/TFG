%
%  CARACTERIZACIÓN DE LA SEÑAL CHIRP R3α
%  TFG: Plataforma de simulación acústica de madera por ultrasonidos
%
%  Reproduce los resultados del Cap. 3 de la tesis (pág. 35-52):
%    - Forma temporal del chirp con ventana gaussiana (Fig.23, pág.43)
%    - Espectro de amplitud y ancho de banda (Fig.24, pág.44)
%    - Comparación con rangos Fco del sensor R3α (Tabla 5, pág.38)
%    - Justificación de la elección del chirp vs. otras señales

clear; clc; close all;

%% Parámetros del chirp R3α (Tabla 5, pág.38)

DT    = 2.0e-7;    % s (mismo que el Par_file)
Ts    = 45e-6;     % s (duración del chirp)
f0_ch = 22e3;      % Hz (frecuencia inicial)
f1_ch = 50e3;      % Hz (frecuencia final)
Fc    = 36e3;      % Hz (frecuencia central = frecuencia de resonancia R3α)
dF    = f1_ch - f0_ch;   % Hz (ancho de barrido)

% Comparación con sensor R6α (Tabla 5)
Ts_R6 = 27e-6;     % s
f0_R6 = 36e3;      % Hz
f1_R6 = 84e3;      % Hz
Fc_R6 = 60e3;      % Hz

fprintf('\n');
fprintf('  CARACTERIZACIÓN SEÑAL CHIRP (Cap. 3 tesis, pág.35-52)\n');
fprintf('\n');

%% Generar chirp R3α

t_ch  = 0:DT:Ts-DT;
mu    = Ts/2;
sigma = Ts/6;
f_t   = (f1_ch-f0_ch)*t_ch/Ts + f0_ch;        % frecuencia instantánea
phase = 2*pi*f_t.*t_ch;                          % fase acumulada
win   = exp(-((t_ch-mu).^2)/(2*sigma^2));        % ventana gaussiana
chirp_R3 = cos(phase) .* win;                    % señal chirp

% Generar chirp R6α (para comparación)
t_R6  = 0:DT:Ts_R6-DT;
f_R6  = (f1_R6-f0_R6)*t_R6/Ts_R6 + f0_R6;
phase_R6 = 2*pi*f_R6.*t_R6;
win_R6   = exp(-((t_R6-Ts_R6/2).^2)/(2*(Ts_R6/6)^2));
chirp_R6 = cos(phase_R6) .* win_R6;

%% Espectros

NSTEP = 5000;
s_R3 = [chirp_R3, zeros(1, NSTEP-length(t_ch))];
s_R6 = [chirp_R6, zeros(1, NSTEP-length(t_R6))];

Nfft  = 2^nextpow2(NSTEP);
f_vec = (0:Nfft/2-1)/(Nfft*DT);

S_R3 = 2*abs(fft(s_R3,Nfft)/NSTEP);  S_R3 = S_R3(1:Nfft/2);
S_R6 = 2*abs(fft(s_R6,Nfft)/NSTEP);  S_R6 = S_R6(1:Nfft/2);

% Ancho de banda -3 dB de cada chirp
thr_R3 = max(S_R3)*10^(-3/20);
thr_R6 = max(S_R6)*10^(-3/20);
BW_R3 = f_vec(find(S_R3>thr_R3,1,'last')) - f_vec(find(S_R3>thr_R3,1,'first'));
BW_R6 = f_vec(find(S_R6>thr_R6,1,'last')) - f_vec(find(S_R6>thr_R6,1,'first'));
Fco_R3_low  = f_vec(find(S_R3>thr_R3,1,'first'))/1e3;
Fco_R3_high = f_vec(find(S_R3>thr_R3,1,'last' ))/1e3;

fprintf('  Sensor R3α (usado en simulación):\n');
fprintf('    f0=%dkHz  f1=%dkHz  Fc=%dkHz  Ts=%dµs\n',f0_ch/1e3,f1_ch/1e3,Fc/1e3,Ts*1e6);
fprintf('    BW -3dB medido:  [%.1f - %.1f] kHz  (BW=%.1f kHz)\n',Fco_R3_low,Fco_R3_high,BW_R3/1e3);
fprintf('    Tesis Tabla 5:   [32.57 - 40.04] kHz  (BW=7.47 kHz)\n\n');

fprintf('  Sensor R6α (referencia tesis):\n');
fprintf('    f0=%dkHz  f1=%dkHz  Fc=%dkHz  Ts=%dµs\n',f0_R6/1e3,f1_R6/1e3,Fc_R6/1e3,Ts_R6*1e6);
fprintf('    BW -3dB medido:  %.1f kHz\n\n', BW_R6/1e3);

%% Frecuencia instantánea

f_inst_R3 = (f1_ch-f0_ch)/Ts * t_ch + f0_ch;   % Hz

%%  FIGURAS DE CARACTERIZACIÓN

% Fig.A: Señal temporal y frecuencia instantánea (reproduce Fig.23, pág.43)
figure('Name','Chirp R3α — Señal temporal','Position',[50 50 850 420]);

subplot(1,2,1); hold on; grid on;
yyaxis left
plot(t_ch*1e6, chirp_R3, 'b-', 'LineWidth',1.5, 'DisplayName','Chirp R3α');
plot(t_ch*1e6, win,      'k--','LineWidth',1.2, 'DisplayName','Ventana gaussiana');
ylabel('Amplitud normalizada');
ylim([-1.2 1.2]);

yyaxis right
plot(t_ch*1e6, f_inst_R3/1e3, 'r-', 'LineWidth',1.5, 'DisplayName','f(t) instantánea');
ylabel('Frecuencia instantánea [kHz]');
ylim([0 80]);

xlabel('Tiempo [µs]');
title('Chirp R3α — forma temporal (Fig.23, pág.43)','FontSize',10,'FontWeight','bold');
legend('Location','north','FontSize',8);
text(Ts*1e6/2, -1.05, ...
    sprintf('f_0=%.0f kHz → f_1=%.0f kHz  |  T_s=%.0f µs  |  F_c=%.0f kHz', ...
    f0_ch/1e3, f1_ch/1e3, Ts*1e6, Fc/1e3), ...
    'FontSize',8, 'HorizontalAlignment','center', 'Color','k');

subplot(1,2,2); hold on; grid on;
plot(t_ch*1e6, chirp_R3,  'b-', 'LineWidth',1.5, 'DisplayName','R3α (22-50 kHz)');
if length(t_R6)*DT <= Ts
    plot(t_R6*1e6, chirp_R6, 'r--','LineWidth',1.5, 'DisplayName','R6α (36-84 kHz)');
end
xlabel('Tiempo [µs]'); ylabel('Amplitud normalizada');
title('Comparación R3α vs. R6α','FontSize',10,'FontWeight','bold');
legend('Location','northeast','FontSize',8);
ylim([-1.2 1.2]);

% Fig.B: Espectros (reproduce Fig.24, pág.44)
figure('Name','Chirp R3α — Espectro','Position',[100 50 850 420]);

subplot(1,2,1); hold on; grid on;
idx_f = f_vec <= 120e3;
plot(f_vec(idx_f)/1e3, S_R3(idx_f), 'b-', 'LineWidth',2, 'DisplayName','R3α');
plot(f_vec(idx_f)/1e3, S_R6(idx_f), 'r--','LineWidth',1.5,'DisplayName','R6α');
xline(32.57,'g:','LineWidth',1.5,'DisplayName','Fco tesis R3α = [32.57-40.04] kHz');
xline(40.04,'g:','LineWidth',1.5,'HandleVisibility','off');
xline(Fc/1e3,'b:','LineWidth',1,'DisplayName',sprintf('F_c R3α = %.0f kHz',Fc/1e3));
xline(Fc_R6/1e3,'r:','LineWidth',1,'DisplayName',sprintf('F_c R6α = %.0f kHz',Fc_R6/1e3));
xlabel('Frecuencia [kHz]'); ylabel('|FFT| normalizada');
title('Espectros de amplitud (Fig.24, pág.44)','FontSize',10,'FontWeight','bold');
legend('Location','northeast','FontSize',7);
xlim([0 110]);

subplot(1,2,2); hold on; grid on;
S_R3_dB = 20*log10(S_R3/max(S_R3));
S_R6_dB = 20*log10(S_R6/max(S_R6));
plot(f_vec(idx_f)/1e3, S_R3_dB(idx_f),'b-', 'LineWidth',2,  'DisplayName','R3α');
plot(f_vec(idx_f)/1e3, S_R6_dB(idx_f),'r--','LineWidth',1.5,'DisplayName','R6α');
yline(-3,'k--','LineWidth',1.5,'DisplayName','-3 dB');
xline(Fco_R3_low, 'b:','LineWidth',1,'HandleVisibility','off');
xline(Fco_R3_high,'b:','LineWidth',1,'HandleVisibility','off');
xlabel('Frecuencia [kHz]'); ylabel('Amplitud [dB]');
title('Espectros en dB — ancho de banda','FontSize',10,'FontWeight','bold');
legend('Location','southwest','FontSize',7);
xlim([0 110]); ylim([-40 5]);
text(mean([Fco_R3_low Fco_R3_high]), -1.5, ...
    sprintf('BW_{-3dB} R3α: %.0f-%.0f kHz', Fco_R3_low, Fco_R3_high), ...
    'FontSize',8,'HorizontalAlignment','center','Color','b');

% Tabla resumen en consola
fprintf('\n');
fprintf('  TABLA: Parámetros del chirp vs. Tesis (Tabla 5, pág.38)\n');
fprintf('\n');
fprintf('  %-22s  %-12s  %-12s\n','Parámetro','Simulación','Tesis R3α');
fprintf('  %s\n',repmat('-',1,50));
tabla = {
    'Sensor',         'R3α',                 'R3α';
    'f_0 [kHz]',      num2str(f0_ch/1e3),    '22';
    'f_1 [kHz]',      num2str(f1_ch/1e3),    '50';
    'F_c [kHz]',      num2str(Fc/1e3),       '36';
    'ΔF [kHz]',       num2str(dF/1e3),       '28';
    'T_s [µs]',       num2str(Ts*1e6),       '45';
    'Fco_low [kHz]',  sprintf('%.1f',Fco_R3_low),  '32.57';
    'Fco_high [kHz]', sprintf('%.1f',Fco_R3_high), '40.04';
    'Ventana',        'Gaussiana',            'Gaussiana';
};
for i=1:size(tabla,1)
    ok = '';
    if strcmp(tabla{i,2},tabla{i,3}), ok='✓'; end
    fprintf('  %-22s  %-12s  %-12s  %s\n',tabla{i,:},ok);
end

fprintf('  NOTA: Fco_tesis = banda del sensor físico R3α (piezoeléctrico)\n');
fprintf('  BW_sim = ancho de banda del chirp matemático generado\n');
fprintf('  Son magnitudes distintas: el chirp se diseña para excitar\n');
fprintf('  la banda del sensor, no para tener el mismo BW que él.\n');