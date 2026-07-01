%
%  BARRIDO DE TAMAÑOS DE DEFECTO
%  Compara dTOF para huecos centrados de 10x10, 5x5 y 2x2 cm
%  vs. caso sano (sin hueco), siguiendo la metodología de la tesis.
%
%  Referencia principal: Fig. 42 y 63, Cap. 4, pág. 62-79
%  "defects located in center of trunk presented larger TOF variations
%   compared to defects in off-centered positions" (pág. 84)
%
%  Configuración:
%    - Fuente: (0.35, 0.50) m   Receptor: (0.65, 0.50) m
%    - Tronco: 30x30 cm  (elementos 30-60 en X y Z)
%    - Defecto centrado: centro en (0.45, 0.45) m
%    - Tamaños: 10x10, 5x5, 2x2 cm
%    - Chirp R3alpha 22-50 kHz (Tabla 5, pág. 38)

clear; clc; close all;

fprintf('\n');
fprintf('  BARRIDO DE TAMAÑOS DE DEFECTO\n');
fprintf('  Comparación con Fig.42/63 de la tesis\n');
fprintf('\n');

%% 1. CONFIGURACIÓN DE CASOS

DT    = 2.0e-7;   % s
NSTEP = 3000;

% --- Directorios de los 4 casos ---
dirs = {
    'Mesh_interna/Pino_sin_hueco_Malla/OUTPUT_FILES/',    'Sin hueco',   0,    NaN;
    'Mesh_interna/Pino_con_hueco_N_C/OUTPUT_FILES/',   '10x10 cm', 10.0,   10.0;
    'Mesh_interna/Pino_con_hueco_N_C_Mediano/OUTPUT_FILES/',     '5x5 cm',    5.0,    5.0;
    'Mesh_interna/Pino_con_hueco_N_C_Pequeño/OUTPUT_FILES/',     '2x2 cm',    2.0,    2.0;
};
% Columnas: directorio, etiqueta, ancho_cm, alto_cm

N_casos = size(dirs,1);

% --- Parámetros físicos ---
dist_SR = 0.30;     % m  ← fuente en xs=0.35, receptor en xr=0.65
V_rad   = 1525;   % m/s  (θ=0°, Christoffel con tus c_ij)
V_tang  = 1120;   % m/s  (θ=90°)
Vp_air  = 343.0;    % m/s
TOF_ref = dist_SR / V_rad * 1e6;   % = 158.16 µs

% --- Señal chirp de referencia (para xcorr) ---
Ts    = 45e-6; f0_ch = 22e3; f1_ch = 50e3;
t_ch  = 0:DT:Ts-DT;
mu = Ts/2; sigma = Ts/6;
chirp_sig = cos(2*pi*((f1_ch-f0_ch)*t_ch/Ts+f0_ch).*t_ch) .* ...
            exp(-((t_ch-mu).^2)/(2*sigma^2));
signal_full = zeros(1,NSTEP);
signal_full(1:length(chirp_sig)) = chirp_sig;

%% 2. CARGA DE TODOS LOS CASOS

fprintf('--- Cargando sismogramas ---\n');

t_all   = cell(N_casos,1);
amp_all = cell(N_casos,1);
ok_all  = false(N_casos,1);

for k = 1:N_casos
    fname = fullfile(dirs{k,1}, 'AA.S0001.BXX.semd');
    if exist(fname,'file')
        datos = load(fname);
        t_k   = datos(:,1);
        amp_k = datos(:,2);
        % Resetear origen temporal (corrige USER_T0 residual)
        t_k   = t_k - t_k(1);
        t_all{k}   = t_k;
        amp_all{k} = amp_k;
        ok_all(k)  = true;
        fprintf('  [%d] %-12s cargado (%d pts)\n', k, dirs{k,2}, length(t_k));
    else
        fprintf('  [%d] %-12s NO ENCONTRADO: %s\n', k, dirs{k,2}, fname);
    end
end

if ~ok_all(1)
    error('Necesito al menos el caso "sin hueco" como referencia base.');
end
if ~any(ok_all(2:end))
    error('Necesito al menos un caso con hueco para comparar.');
end

%% 3. DETECCIÓN DE TOF EN CADA CASO

fprintf('\n--- Detección de TOF (correlación cruzada, Cap. 3) ---\n');

TOF_xcorr = NaN(N_casos,1);
TOF_thresh = NaN(N_casos,1);
TOF_env    = NaN(N_casos,1);   

for k = 1:N_casos
    if ~ok_all(k), continue; end
    [TOF_xcorr(k), TOF_thresh(k), TOF_env(k)] = detectar_TOF_local(...   % TOF_env
        amp_all{k}, signal_full, DT, NSTEP);
    fprintf('  [%d] %-12s  xcorr=%.2f µs  threshold=%.2f µs  env=%.2f µs\n', ...
        k, dirs{k,2}, TOF_xcorr(k), TOF_thresh(k), TOF_env(k));   % env
end

% Referencia: caso sin hueco medido por envolvente (más preciso para TOF absoluto)
TOF_base = TOF_env(1);   % ← antes era TOF_xcorr(1)
fprintf('\n  TOF base (sin hueco, envolvente): %.2f µs\n', TOF_base);
fprintf('  TOF referencia analítico:         %.2f µs\n', TOF_ref);
fprintf('  Error modelo vs. medición:        %.2f%%\n', ...
    abs(TOF_base - TOF_ref)/TOF_ref*100);

%% 4. CALCULAR dTOF PARA CADA TAMAÑO

fprintf('\n--- dTOF respecto al caso sano (método: envolvente) ---\n');
fprintf('  %-12s  %-10s  %-10s  %-8s  %-12s\n', ...
    'Caso','TOF_env[µs]','dTOF[µs]','dTOF[%]','Ref tesis');
fprintf('  %s\n', repmat('-',1,60));

% Referencia tesis Fig.63 (pág.79) — Oak, sensor radial opuesto:
%   Defecto Ø 2.9 cm:  dTOF_sim= 13µs (+11%), dTOF_exp= 18µs (+12%)
%   Defecto Ø 5.1 cm:  dTOF_sim= 26µs (+22%), dTOF_exp= 50µs (+39%)
%   Defecto Ø 7.6 cm:  dTOF_sim= 46µs (+34%), dTOF_exp= 75µs (+50%)
% Tus tamaños equivalentes aproximados:
%   2x2 cm  ≈ Ø2.9 cm → ref: +11% (sim)
%   5x5 cm  ≈ Ø5.1 cm → ref: +22% (sim)
%   10x10cm ≈ Ø10 cm  → ref: +34% (sim, Cap.5 numérico)

ref_tesis = {'base','~+34% (sim)','~+22% (sim)','~+11% (sim)'};

dTOF_us  = NaN(N_casos,1);
dTOF_pct = NaN(N_casos,1);

for k = 1:N_casos
    if ~ok_all(k), continue; end
    dTOF_us(k)  = TOF_env(k) - TOF_base;
    dTOF_pct(k) = dTOF_us(k) / TOF_base * 100;
    fprintf('  %-12s  %10.2f  %10.2f  %7.1f%%  %s\n', ...
        dirs{k,2}, TOF_env(k), dTOF_us(k), dTOF_pct(k), ref_tesis{k});
end

%% 5. ESTIMACIÓN ANALÍTICA DE dTOF (rayo recto, si pasa por el hueco)

fprintf('\n--- Estimación analítica dTOF (si rayo atravesara el defecto) ---\n');
fprintf('  El rayo (z=0.45) ATRAVIESA el hueco centrado.\n');
fprintf('  La estimación analítica asume que el rayo cruza el hueco en línea recta.\n');
fprintf('  El campo de ondas real (SPECFEM2D) incluye difracción → dTOF real < analítico.\n\n');

xs = 0.35; xr = 0.65;
hueco_cx = 0.50;   % centro x del hueco (fijo)

for k = 2:N_casos
    if isnan(dirs{k,4}), continue; end
    d_cm   = dirs{k,3};
    d_m    = d_cm / 100;
    d_madera = dist_SR - d_m;           % trayectoria fuera del hueco
    TOF_an_sin = dist_SR / V_rad * 1e6; % µs
    TOF_an_con = (d_madera/V_rad + d_m/Vp_air) * 1e6;
    dTOF_an = TOF_an_con - TOF_an_sin;
    dTOF_an_pct = dTOF_an / TOF_an_sin * 100;
    fprintf('  %-12s: TOF_sin=%.2f µs, TOF_con=%.2f µs, dTOF=%.2f µs (+%.1f%%)\n',...
        dirs{k,2}, TOF_an_sin, TOF_an_con, dTOF_an, dTOF_an_pct);
end

%% 6. TABLA RESUMEN DE VALIDACIÓN

fprintf('\n');
fprintf('  TABLA RESUMEN BARRIDO DE TAMAÑOS\n');
fprintf('\n');
fprintf('  %-12s  %-10s  %-10s  %-12s  %-14s\n',...
    'Caso','TOF[µs]','dTOF[µs]','dTOF[%]','Ref. tesis[%]');
fprintf('  %s\n', repmat('-',1,64));
refs_pct = {'base', '~34', '~22', '~11'};
for k = 1:N_casos
    if ~ok_all(k), continue; end
    fprintf('  %-12s  %10.2f  %10.2f  %11.1f%%  %s\n',...
        dirs{k,2}, TOF_env(k), dTOF_us(k), dTOF_pct(k), refs_pct{k});
end

fprintf('\n  Tendencia esperada (tesis Fig.42, pág.62):\n');
fprintf('  → dTOF aumenta con el tamaño del defecto ✓ si se cumple\n');
fprintf('  → El defecto mayor (10x10) debe dar el dTOF más alto\n');
fprintf('  → Si tu dTOF es pequeño (<<+10%%), el rayo es tangente al hueco:\n');
fprintf('    mueve el hueco 5 cm hacia abajo (zmin: 30→25) para que el\n');
fprintf('    rayo lo atraviese y obtengas efectos comparables a Fig.63.\n');

%% 7. FIGURAS

colores = {'k','b','r','m'};
estilos = {'-','--',':','-.'};

% Figura 1: Sismogramas superpuestos
figure('Name','Tamaños: Sismogramas','Position',[50 50 950 550]);
hold on; grid on;
for k = 1:N_casos
    if ~ok_all(k), continue; end
    plot(t_all{k}*1e6, amp_all{k}/max(abs(amp_all{k})), ...
        [colores{k} estilos{k}], 'LineWidth',1.5, 'DisplayName', dirs{k,2});
end
if ~isnan(TOF_base)
    xline(TOF_base,'k--','LineWidth',1,'DisplayName',sprintf('TOF base=%.1f µs',TOF_base));
end
xline(TOF_ref,'g-.','LineWidth',1,'DisplayName',sprintf('TOF analítico=%.1f µs',TOF_ref));
xlabel('Tiempo [µs]'); ylabel('Amplitud normalizada');
title('Sismogramas: pino con huecos de distintos tamaños','FontSize',11);
legend('Location','northeast','FontSize',8);
Tlim = min(t_all{1}(end)*1e6, 800);
xlim([0 Tlim]);

% Figura 2: dTOF vs. tamaño del defecto
idx_validos = find(ok_all & ~isnan(dTOF_us) & (1:N_casos)' > 1);
if length(idx_validos) >= 2
    tamanios_cm = zeros(1, length(idx_validos));
    dTOF_plot   = zeros(1, length(idx_validos));
    for j = 1:length(idx_validos)
        k = idx_validos(j);
        tamanios_cm(j) = dirs{k,3};
        dTOF_plot(j)   = dTOF_pct(k);
    end

    figure('Name','Tamaños: dTOF vs. tamaño','Position',[100 100 700 450]);
    subplot(1,2,1);
    bar(tamanios_cm, dTOF_us(idx_validos), 0.5, 'FaceColor',[0.3 0.6 0.9]);
    hold on; grid on;
    % Barras de referencia tesis (sim)
    ref_us_tesis = [NaN, 13, 26, 46];   % tamaños 0,2,5,10 → NaN para "sin hueco"
    scatter([2 5 10], [13 26 46], 80, 'r^', 'filled', 'DisplayName','Tesis sim (Fig.63)');
    xlabel('Tamaño defecto [cm]'); ylabel('dTOF [µs]');
    title('dTOF vs. tamaño (ref. Fig.63, pág.79)','FontSize',10);
    legend({'Tu simulación','Tesis (sim)'},'Location','northwest','FontSize',8);

    subplot(1,2,2);
    bar(tamanios_cm, dTOF_plot, 0.5, 'FaceColor',[0.3 0.8 0.5]);
    hold on; grid on;
    scatter([2 5 10], [11 22 34], 80, 'r^', 'filled', 'DisplayName','Tesis sim (%)');
    scatter([2 5 10], [12 39 50], 80, 'bs', 'filled', 'DisplayName','Tesis exp (%)');
    xlabel('Tamaño defecto [cm]'); ylabel('dTOF [%]');
    title('dTOF [%] vs. tamaño','FontSize',10);
    legend('Location','northwest','FontSize',8);
    yline(0,'k-','LineWidth',0.5);
end

% Figura 3: Correlación entre casos
figure('Name','Tamaños: Correlaciones cruzadas','Position',[150 150 900 380]);
hold on; grid on;
for k = 2:N_casos
    if ~ok_all(k) || ~ok_all(1), continue; end
    s_base_n = amp_all{1} / max(abs(amp_all{1}));
    s_k_n    = amp_all{k} / max(abs(amp_all{k}));
    N_min    = min(length(s_base_n), length(s_k_n));
    [cc, ll] = xcorr(s_base_n(1:N_min), s_k_n(1:N_min), 'coeff');
    plot(ll*DT*1e6, cc, [colores{k} estilos{k}], 'LineWidth',1.5, ...
        'DisplayName', sprintf('%s vs sin hueco', dirs{k,2}));
end
xlabel('Desfase [µs]'); ylabel('Correlación cruzada normalizada');
title('Correlación de cada caso con defecto vs. caso sano','FontSize',10);
legend('Location','northeast','FontSize',8);
xlim([-200 200]);

%% FUNCIONES LOCALES

function [TOF_xc_us, TOF_th_us, TOF_en_us] = detectar_TOF_local(amp, signal_full, DT, NSTEP)
% Detecta TOF por tres métodos con ventana física [130µs, 500µs]
% La ventana evita que la xcorr elija picos espurios fuera del rango físico

    T_total = NSTEP * DT;
    N       = length(amp);
    TOF_min_s = 0.130e-3;               % 130 µs mínimo físico
    TOF_max_s = min(0.500e-3, T_total); % 500 µs máximo razonable

    % 1. xcorr — solo el chirp (sin ceros de relleno), ventana física
    idx_nz = find(signal_full ~= 0);
    ch_n   = signal_full(idx_nz(1):idx_nz(end))';
    ch_n   = ch_n / (norm(ch_n) + eps);
    rec_n  = amp(:) / (norm(amp(:)) + eps);
    [xc, lags] = xcorr(rec_n, ch_n, round(T_total/DT));
    lags_s = lags * DT;
    v = (lags_s >= TOF_min_s) & (lags_s <= TOF_max_s);
    if ~any(v), v = (lags_s >= 0) & (lags_s <= T_total); end
    lv = lags_s(v); xv = xc(v);
    [~, idx] = max(xv);
    TOF_xc_us = lv(idx) * 1e6;

    % 2. Threshold 5% del máximo, ventana física
    i0  = round(TOF_min_s / DT);
    i1  = min(round(TOF_max_s / DT), N);
    seg = amp(i0:i1);
    ith = find(abs(seg) > 0.05*max(abs(amp)), 1, 'first');
    if ~isempty(ith)
        TOF_th_us = (ith + i0 - 1) * DT * 1e6;
    else
        TOF_th_us = NaN;
    end

    % 3. Envolvente RMS local (más fiable con hueco grande)
    %    No requiere Signal Processing Toolbox
    vw  = max(round(1/(36e3*DT)), 10);   % ventana ≈ 1 período de Fc=36kHz
    env = sqrt(movmean(amp(:).^2, vw));
    seg_e = env(i0:i1);
    ie = find(seg_e > 0.02*max(env), 1, 'first');
    if ~isempty(ie)
        TOF_en_us = (ie + i0 - 1) * DT * 1e6;
    else
        TOF_en_us = NaN;
    end
end
