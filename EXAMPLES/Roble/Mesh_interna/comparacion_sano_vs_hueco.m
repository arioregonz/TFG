% 
%  COMPARACIÓN SANO vs. HUECO — Métricas cuantitativas
%  TFG: Plataforma de simulación acústica de madera por ultrasonidos
%
%  Calcula y compara:
%    - TOF por umbral y envolvente (con ventana física)
%    - dTOF entre casos (método xcorr entre señales)
%    - Error RMS, energía y amplitud
%    - Correlación temporal alineada y correlación espectral
%
%  Complementa validacion_par_file.m añadiendo métricas de similitud
%  entre señales que permiten cuantificar el efecto del defecto.

clear; clc; close all;

fprintf('\n');
fprintf(' COMPARACION SANO vs HUECO\n');
fprintf('\n');

%% 1. PARÁMETROS Y RUTAS
DT      = 2.0e-7;     % s
TOF_ref = 158.16e-6;  % s  (modelo analítico: 0.30m / 1896.8 m/s)

file_sano  = 'Mesh_interna/Roble_sin_hueco_Norm/OUTPUT_FILES/AA.S0001.BXX.semd';
file_hueco = 'Mesh_interna/Roble_con_hueco_N_C/OUTPUT_FILES/AA.S0001.BXX.semd';

%% 2. CARGAR Y NORMALIZAR
A = load(file_sano);
B = load(file_hueco);

% Resetear origen temporal (corrige USER_T0 residual)
t1 = A(:,1) - A(1,1);
t2 = B(:,1) - B(1,1);
s1_raw = A(:,2);
s2_raw = B(:,2);

% Normalizar por el máximo de amplitud
s1 = s1_raw / max(abs(s1_raw));
s2 = s2_raw / max(abs(s2_raw));

fprintf('Archivos cargados:\n');
fprintf('  Sano:  %d puntos, T_total=%.0f µs\n', length(t1), t1(end)*1e6);
fprintf('  Hueco: %d puntos, T_total=%.0f µs\n\n', length(t2), t2(end)*1e6);

%% 3. DETECCIÓN DE TOF CON VENTANA FÍSICA [130µs, 500µs]
%  La ventana evita que el umbral se dispare por ruido numérico inicial
%  y que la xcorr elija picos espurios fuera del rango físico posible.

TOF_min_s = 0.130e-3;   % 130 µs  (ninguna onda puede llegar antes)
TOF_max_s = 0.500e-3;   % 500 µs  (límite razonable para 0.30m de madera)

i0 = round(TOF_min_s / DT);
i1 = min(round(TOF_max_s / DT), length(s1));

% Umbral 5% del máximo (más robusto que umbral fijo 0.15)
umbral = 0.05;
idx1 = find(abs(s1(i0:i1)) > umbral, 1, 'first');
idx2 = find(abs(s2(i0:i1)) > umbral, 1, 'first');

TOF_sano  = (idx1 + i0 - 1) * DT;
TOF_hueco = (idx2 + i0 - 1) * DT;
dTOF_thr  = TOF_hueco - TOF_sano;

% Envolvente RMS (más fiable para hueco grande)
vw   = max(round(1/(36e3*DT)), 10);
env1 = sqrt(movmean(s1.^2, vw));
env2 = sqrt(movmean(s2.^2, vw));
ie1  = find(env1(i0:i1) > 0.02*max(env1), 1, 'first');
ie2  = find(env2(i0:i1) > 0.02*max(env2), 1, 'first');
TOF_env_sano  = (ie1 + i0 - 1) * DT;
TOF_env_hueco = (ie2 + i0 - 1) * DT;
dTOF_env      = TOF_env_hueco - TOF_env_sano;

fprintf('--- Detección de TOF ---\n');
fprintf('  %-12s  %-14s  %-14s  %s\n', 'Método', 'TOF sano [µs]', 'TOF hueco [µs]', 'dTOF [µs] (%)');
fprintf('  %s\n', repmat('-',1,60));
fprintf('  %-12s  %-14.2f  %-14.2f  %.2f (%.1f%%)\n', 'Umbral 5%%', ...
    TOF_sano*1e6, TOF_hueco*1e6, dTOF_thr*1e6, dTOF_thr/TOF_sano*100);
fprintf('  %-12s  %-14.2f  %-14.2f  %.2f (%.1f%%)\n', 'Envolvente', ...
    TOF_env_sano*1e6, TOF_env_hueco*1e6, dTOF_env*1e6, dTOF_env/TOF_env_sano*100);
fprintf('  TOF referencia analítico: %.2f µs\n', TOF_ref*1e6);
fprintf('  Error envolvente vs ref:  %.1f%%\n\n', ...
    abs(TOF_env_sano - TOF_ref)/TOF_ref*100);

%% 4. dTOF POR XCORR ENTRE SEÑALES RECIBIDAS
%  xcorr(s2,s1) mide el retardo de s2 respecto a s1 → dTOF directo
%  Usar ventana física para evitar picos espurios

N = min(length(s1), length(s2));
[c, lags] = xcorr(s2(1:N), s1(1:N), 'coeff');
lags_s     = lags * DT;

% Buscar solo en ventana de dTOF físicamente posible: [-200µs, +400µs]
v_xcorr = (lags_s >= -0.200e-3) & (lags_s <= 0.400e-3);
[val_xc, idx_xc] = max(c(v_xcorr));
lags_v    = lags_s(v_xcorr);
dTOF_xcorr = lags_v(idx_xc);

fprintf('--- dTOF por xcorr entre señales recibidas ---\n');
fprintf('  Pico de correlación: %.4f\n', val_xc);
fprintf('  dTOF_xcorr = %.2f µs\n\n', dTOF_xcorr*1e6);

%% 5. MÉTRICAS DE SIMILITUD

% Error RMS (sobre señales normalizadas, ventana física [130µs, 500µs])
s1_win = s1(i0:i1);
s2_win = s2(i0:min(i1, length(s2)));
N_win  = min(length(s1_win), length(s2_win));
error_rms = sqrt(mean((s1_win(1:N_win) - s2_win(1:N_win)).^2));

% Energía (señales sin normalizar, trapezoidal)
E1 = trapz(t1, s1_raw.^2);
E2 = trapz(t2, s2_raw.^2);
ratio_E = E2 / E1;

% Amplitud máxima
A1 = max(abs(s1_raw));
A2 = max(abs(s2_raw));
ratio_A = A2 / A1;

% Correlación temporal con alineación por desplazamiento (no circshift)
% Alinear s2 respecto a s1 usando el lag de la xcorr
lag_muestras = round(dTOF_xcorr / DT);
if lag_muestras >= 0
    s2_alin = [zeros(lag_muestras,1); s2(1:N-lag_muestras)];
else
    s2_alin = [s2(1-lag_muestras:N); zeros(-lag_muestras,1)];
end
R_temp = corrcoef(s1(1:N), s2_alin(1:N));
corr_alineada = R_temp(1,2);

% Correlación espectral (solo en la banda útil del chirp: 15-80 kHz)
Nfft   = 2^nextpow2(N);
f_vec  = (0:Nfft/2-1) / (Nfft * DT);
S1_f   = abs(fft(s1(1:N), Nfft));
S2_f   = abs(fft(s2(1:N), Nfft));
S1_f   = S1_f(1:Nfft/2);
S2_f   = S2_f(1:Nfft/2);

% Banda útil: 15-80 kHz
idx_band = (f_vec >= 15e3) & (f_vec <= 80e3);
R_spec  = corrcoef(S1_f(idx_band), S2_f(idx_band));
corr_espectral = R_spec(1,2);

fprintf('--- Métricas de similitud ---\n');
fprintf('  Error RMS (ventana física):   %.4f\n', error_rms);
fprintf('  Ratio de energía E_hueco/E_sano: %.4f  (%.1f%%)\n', ratio_E, ratio_E*100);
fprintf('  Ratio amplitud A_hueco/A_sano:   %.4f  (%.1f%%)\n', ratio_A, ratio_A*100);
fprintf('  Correlación temporal alineada:   %.4f\n', corr_alineada);
fprintf('  Correlación espectral (15-80kHz):%.4f\n\n', corr_espectral);

%% 6. RESUMEN Y DIAGNÓSTICO

fprintf('\n');
fprintf('RESUMEN\n');
fprintf('\n');
fprintf('  dTOF umbral       = %+.2f µs (%+.1f%%)\n', dTOF_thr*1e6, dTOF_thr/TOF_sano*100);
fprintf('  dTOF envolvente   = %+.2f µs (%+.1f%%)  ← más fiable\n', dTOF_env*1e6, dTOF_env/TOF_env_sano*100);
fprintf('  dTOF xcorr        = %+.2f µs\n', dTOF_xcorr*1e6);
fprintf('  Corr. alineada    = %.4f\n', corr_alineada);
fprintf('  Corr. espectral   = %.4f\n', corr_espectral);
fprintf('  Error RMS         = %.4f\n', error_rms);
fprintf('  Ratio energía     = %.4f\n', ratio_E);
fprintf('  Ratio amplitud    = %.4f\n', ratio_A);

% Diagnóstico automático
fprintf('\n  DIAGNÓSTICO:\n');

if dTOF_env*1e6 > 10
    fprintf('  → dTOF significativo (>10µs): defecto detectable ✓\n');
else
    fprintf('  → dTOF pequeño (<10µs): defecto marginal o rayo no lo cruza\n');
end

if corr_alineada > 0.90
    fprintf('  → Alta correlación temporal (%.2f): señal con hueco mantiene forma de chirp\n', corr_alineada);
elseif corr_alineada > 0.70
    fprintf('  → Correlación media (%.2f): hueco distorsiona parcialmente la señal\n', corr_alineada);
else
    fprintf('  → Baja correlación (%.2f): hueco bloquea camino directo\n', corr_alineada);
    fprintf('     La señal recibida son ondas difractadas (forma ≠ chirp emitido)\n');
end

if ratio_E < 0.5
    fprintf('  → Energía hueco << sano (%.0f%%): el defecto atenúa fuertemente\n', ratio_E*100);
elseif ratio_E < 0.80
    fprintf('  → Energía hueco < sano (%.0f%%): atenuación moderada por el defecto\n', ratio_E*100);
else
    fprintf('  → Energía similar (%.0f%%): el defecto no atenúa mucho\n', ratio_E*100);
end

%% 7. FIGURAS

% Fig.1: Sismogramas superpuestos con marcadores de TOF
figure('Name','Sano vs Hueco — Sismogramas','Position',[50 50 900 500]);

subplot(2,1,1); hold on; grid on;
plot(t1*1e6, s1, 'b-', 'LineWidth',1.5, 'DisplayName','Sano');
plot(t2*1e6, s2, 'r-', 'LineWidth',1.5, 'DisplayName','Con hueco');
xline(TOF_ref*1e6,        'g--', 'LineWidth',1.5, 'DisplayName', sprintf('TOF_{ref}=%.0fµs',TOF_ref*1e6));
xline(TOF_env_sano*1e6,   'b:',  'LineWidth',2.0, 'DisplayName', sprintf('TOF_{env,sano}=%.0fµs',TOF_env_sano*1e6));
xline(TOF_env_hueco*1e6,  'r:',  'LineWidth',2.0, 'DisplayName', sprintf('TOF_{env,hueco}=%.0fµs',TOF_env_hueco*1e6));
xlabel('Tiempo [µs]'); ylabel('Amplitud norm.');
title('Sismogramas normalizados — sano vs. con hueco','FontSize',10,'FontWeight','bold');
legend('Location','northeast','FontSize',7);
xlim([0 t1(end)*1e6]);

subplot(2,1,2); hold on; grid on;
plot(t1*1e6, s1, 'b-', 'LineWidth',1.8, 'DisplayName','Sano');
plot(t2*1e6, s2, 'r-', 'LineWidth',1.8, 'DisplayName','Con hueco');
xline(TOF_ref*1e6,       'g--', 'LineWidth',1.5);
xline(TOF_env_sano*1e6,  'b:',  'LineWidth',2.0);
xline(TOF_env_hueco*1e6, 'r:',  'LineWidth',2.0);
xlabel('Tiempo [µs]'); ylabel('Amplitud norm.');
title(sprintf('ZOOM primera llegada | dTOF_{env} = %+.0f µs (%+.1f%%)', ...
    dTOF_env*1e6, dTOF_env/TOF_env_sano*100), 'FontSize',10,'FontWeight','bold');
legend('Location','northeast','FontSize',7);
xlim([max(0, TOF_ref*1e6 - 50), TOF_ref*1e6 + 300]);

% Fig.2: xcorr entre señales recibidas
figure('Name','Correlación cruzada entre señales','Position',[100 100 800 380]);
hold on; grid on;
plot(lags_s*1e6, c, 'k-', 'LineWidth',1.5);
xline(dTOF_xcorr*1e6, 'r-', 'LineWidth',2, ...
    'DisplayName', sprintf('dTOF_{xcorr}=%.0f µs', dTOF_xcorr*1e6));
xline(0, 'g--', 'LineWidth',1, 'DisplayName', 'Retardo cero');
xlabel('Retardo [µs]'); ylabel('Correlación normalizada');
title(sprintf('xcorr(hueco, sano) — dTOF = %.0f µs', dTOF_xcorr*1e6), ...
    'FontSize',10,'FontWeight','bold');
legend('FontSize',8);
xlim([-200 500]);

% Fig.3: Espectros en banda útil
figure('Name','Comparación espectral','Position',[150 150 800 380]);

subplot(1,2,1); hold on; grid on;
plot(f_vec/1e3, S1_f/max(S1_f), 'b-', 'LineWidth',1.5, 'DisplayName','Sano');
plot(f_vec/1e3, S2_f/max(S2_f), 'r-', 'LineWidth',1.5, 'DisplayName','Con hueco');
xline(22, 'k:', 'LineWidth',1, 'DisplayName','f_0=22 kHz');
xline(50, 'k:', 'LineWidth',1, 'HandleVisibility','off');
xlabel('Frecuencia [kHz]'); ylabel('Magnitud norm.');
title('Espectros — banda del chirp','FontSize',10,'FontWeight','bold');
legend('Location','northeast','FontSize',8);
xlim([0 80]);

subplot(1,2,2); hold on; grid on;
plot(f_vec(idx_band)/1e3, S1_f(idx_band)/max(S1_f(idx_band)), 'b-', 'LineWidth',1.5);
plot(f_vec(idx_band)/1e3, S2_f(idx_band)/max(S2_f(idx_band)), 'r-', 'LineWidth',1.5);
xlabel('Frecuencia [kHz]'); ylabel('Magnitud norm.');
title(sprintf('Banda 15-80 kHz | corr. esp. = %.4f', corr_espectral), ...
    'FontSize',10,'FontWeight','bold');
legend({'Sano','Con hueco'},'FontSize',8);
xlim([15 80]);

%% 8. GUARDAR
save('comparacion_sano_hueco.mat', ...
    'TOF_env_sano', 'TOF_env_hueco', 'dTOF_env', ...
    'TOF_sano', 'TOF_hueco', 'dTOF_thr', ...
    'dTOF_xcorr', 'val_xc', ...
    'error_rms', 'ratio_E', 'ratio_A', ...
    'corr_alineada', 'corr_espectral', ...
    't1', 's1', 't2', 's2');
fprintf('\n  Resultados guardados: comparacion_sano_hueco.mat\n');