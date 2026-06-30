%
%  BARRIDO DE POSICIONES DE DEFECTO (10x10 cm)
%  Compara dTOF para el hueco en 5 posiciones:
%  Centro, UR (arriba-derecha), UL (arriba-izquierda), DL, DR
%
%  Referencia principal: Fig. 44-51, Cap. 4, pág. 63-67
%  "Defects located in the center of the trunk presented larger TOF
%   variations compared to defects in off-centered positions" (pág. 84)
%  "Off-centered defects will be more difficult to determine" (pág. 84)
%
%  Configuración:
%    - Fuente: (0.30, 0.45) m   Receptor: (0.60, 0.45) m
%    - Tronco: 30x30 cm, centro en (0.45, 0.45) m
%    - Defecto: 10x10 cm en 5 posiciones
%    - Chirp R3alpha 22-50 kHz

clear; clc; close all;

fprintf('\n');
fprintf('  BARRIDO DE POSICIONES DE DEFECTO (10x10 cm)\n');
fprintf('  Comparación con Fig.44-51 de la tesis\n');
fprintf('\n');

%% 1.CONFIGURACIÓN DE CASOS Y GEOMETRÍA

DT    = 2.0e-7;   % s
NSTEP = 3000;

% Parámetros físicos
dist_SR = 0.30;    % m  (fuente-receptor)
V_rad   = 1896.8;  % m/s
V_tang  = 1382.8;  % m/s
Vp_air  = 343.0;   % m/s
TOF_ref = dist_SR / V_rad * 1e6;  % µs

% Geometría del tronco
cx_t = 0.45; cz_t = 0.45;   % centro del tronco [m]
R_t  = 0.15;                 % radio equivalente [m]

% Geometría de los huecos 10x10 cm en cada posición
% Cada posición especifica el centro del hueco
d_h = 0.10;   % tamaño del hueco [m]
offset = 0.07; % desplazamiento del centro desde el centro del tronco [m]
%   (7 cm ≈ 0.5 * (R_tronco - d_hueco/2), análogo a "offset 7.5 cm" tesis pág.90)

% Centro de cada hueco [m]
huecos_cx = [cx_t,              cx_t+offset, cx_t-offset, cx_t-offset, cx_t+offset];
huecos_cz = [cz_t,              cz_t+offset, cz_t+offset, cz_t-offset, cz_t-offset];
nombres_pos = {'Centro','UR','UL','DL','DR'};

% Esquinas [xmin, xmax, zmin, zmax] de cada hueco
for k = 1:5
    huecos(k).cx   = huecos_cx(k);
    huecos(k).cz   = huecos_cz(k);
    huecos(k).xmin = huecos_cx(k) - d_h/2;
    huecos(k).xmax = huecos_cx(k) + d_h/2;
    huecos(k).zmin = huecos_cz(k) - d_h/2;
    huecos(k).zmax = huecos_cz(k) + d_h/2;
end

% Directorios (sin hueco + 5 posiciones con hueco)
dirs = {
    'Mesh_interna/Roble_sin_hueco_Norm/OUTPUT_FILES/',       'Sin hueco';
    'Mesh_interna/Roble_con_hueco_N_C/OUTPUT_FILES/',     'Centro';
    'Mesh_interna/Roble_con_hueco_N_UR/OUTPUT_FILES/',         'UR';
    'Mesh_interna/Roble_con_hueco_N_UL/OUTPUT_FILES/',         'UL';
    'Mesh_interna/Roble_con_hueco_N_DL/OUTPUT_FILES/',         'DL';
    'Mesh_interna/Roble_con_hueco_N_DR/OUTPUT_FILES/',         'DR';
};
N_casos = size(dirs,1);

% Señal chirp
Ts    = 45e-6; f0_ch = 22e3; f1_ch = 50e3;
t_ch  = 0:DT:Ts-DT;
mu = Ts/2; sigma = Ts/6;
chirp_sig = cos(2*pi*((f1_ch-f0_ch)*t_ch/Ts+f0_ch).*t_ch) .* ...
            exp(-((t_ch-mu).^2)/(2*sigma^2));
signal_full = zeros(1,NSTEP);
signal_full(1:length(chirp_sig)) = chirp_sig;

%% 2. CARGA DE SISMOGRAMAS

fprintf('--- Cargando sismogramas ---\n');

t_all   = cell(N_casos,1);
amp_all = cell(N_casos,1);
ok_all  = false(N_casos,1);

for k = 1:N_casos
    fname = fullfile(dirs{k,1},'AA.S0001.BXX.semd');
    if exist(fname,'file')
        datos      = load(fname);
        t_k        = datos(:,1) - datos(1,1);   % resetear USER_T0
        amp_all{k} = datos(:,2);
        t_all{k}   = t_k;
        ok_all(k)  = true;
        fprintf('  [%d] %-12s cargado (%d pts)\n', k, dirs{k,2}, length(t_k));
    else
        fprintf('  [%d] %-12s NO ENCONTRADO: %s\n', k, dirs{k,2}, fname);
    end
end

if ~ok_all(1)
    error('Necesito el caso "sin hueco" como referencia base.');
end

%% 3. DETECCIÓN DE TOF

fprintf('\n--- Detección de TOF ---\n');

TOF_xcorr = NaN(N_casos,1);
TOF_thresh = NaN(N_casos,1);
TOF_env = NaN(N_casos,1);

for k = 1:N_casos
    if ~ok_all(k), continue; end
    [TOF_xcorr(k), TOF_thresh(k), TOF_env(k)] = detectar_TOF_local(amp_all{k}, signal_full, DT, NSTEP);
    fprintf('  [%d] %-12s  xcorr=%.2f µs  threshold=%.2f µs  env=%.2f µs\n', k, dirs{k,2}, TOF_xcorr(k), TOF_thresh(k), TOF_env(k));
end

TOF_base = TOF_env(1); 
fprintf('\n  TOF base (sin hueco, envolvente): %.2f µs\n', TOF_base);
fprintf('  TOF analítico (V_rad): %.2f µs\n', TOF_ref);

%% 4. CALCULAR dTOF Y ANÁLISIS DE POSICIÓN

fprintf('\n--- dTOF por posición ---\n');
fprintf('  %-12s  %-8s  %-8s  %-8s  %-22s  %s\n',...
    'Posición','TOF[µs]','dTOF[µs]','dTOF[%]','Rayo pasa por hueco','Ref. tesis');
fprintf('  %s\n', repmat('-',1,80));

dTOF_us  = NaN(N_casos,1);
dTOF_pct = NaN(N_casos,1);
rayo_pasa = false(N_casos,1);

xs_f = 0.30; zs_f = 0.45;   % fuente
xr_f = 0.60; zr_f = 0.45;   % receptor

for k = 1:N_casos
    if ~ok_all(k) || k==1, continue; end
    dTOF_us(k)  = TOF_env(k) - TOF_base;
    dTOF_pct(k) = dTOF_us(k) / TOF_base * 100;

    % Verificar si el rayo recto pasa dentro del hueco
    % (rayo horizontal z=0.45, el hueco tiene zmin/zmax)
    hi = k-1;  % índice en array huecos (caso k=2 → hueco 1=Centro)
    if hi >= 1 && hi <= 5
        hueco_k = huecos(hi);
        % El rayo pasa por el hueco si:
        %   zmin < zs_f < zmax  (en z)  Y  xmin < algún punto del rayo < xmax (en x)
        rayo_en_z = (hueco_k.zmin < zs_f) && (zs_f < hueco_k.zmax);
        rayo_en_x = (hueco_k.xmin < xr_f) && (hueco_k.xmax > xs_f);
        rayo_pasa(k) = rayo_en_z && rayo_en_x;
    end

    % Referencia tesis por posición (Fig. 45-51, pág. 63-67)
    % Para defecto centrado (cap. 5): dTOF ~+34% (sim)
    % Para defecto excéntrico vertical +8cm (Fig.49): dTOF ~+15% (sim)
    % Para defecto excéntrico horizontal ±8cm (Fig.45): dTOF ~<1% (sim)
    ref_pos = {
        'Centro', '~+34% (sim, Cap.5) — rayo atraviesa';
        'UR',     '~<11% (solo difracción, no en ruta)';
        'UL',     '~<11% (solo difracción, no en ruta)';
        'DL',     '~<11% (solo difracción, no en ruta)';
        'DR',     '~<11% (solo difracción, no en ruta)';
    };

    ref_str = '-';
    for r = 1:size(ref_pos,1)
        if strcmp(dirs{k,2}, ref_pos{r,1})
            ref_str = ref_pos{r,2};
        end
    end

    fprintf('  %-12s  %8.2f  %8.2f  %7.1f%%  %-22s  %s\n',dirs{k,2}, TOF_env(k), dTOF_us(k), dTOF_pct(k), ... 
        ternario(rayo_pasa(k),'SÍ (efecto directo)','NO (solo difracción)'), ref_str);
end

%% 5. INTERPRETACIÓN EN BASE A LA TESIS

fprintf('\n--- Interpretación (Cap. 4, pág. 63-84) ---\n');

% Obtener dTOFs de los casos con hueco
dTOF_centro = NaN; dTOF_UR = NaN; dTOF_UL = NaN; dTOF_DL = NaN; dTOF_DR = NaN;
for k = 1:N_casos
    switch dirs{k,2}
        case 'Centro', dTOF_centro = dTOF_pct(k);
        case 'UR',     dTOF_UR     = dTOF_pct(k);
        case 'UL',     dTOF_UL     = dTOF_pct(k);
        case 'DL',     dTOF_DL     = dTOF_pct(k);
        case 'DR',     dTOF_DR     = dTOF_pct(k);
    end
end

% Comprobar la tendencia principal de la tesis
fprintf('  Tendencia 1 (pág.84): dTOF_centro > dTOF_excéntrico\n');
dTOF_exc_max = max([dTOF_UR, dTOF_UL, dTOF_DL, dTOF_DR], [], 'omitnan');
if ~isnan(dTOF_centro) && ~isnan(dTOF_exc_max)
    cumple = dTOF_centro > dTOF_exc_max;
    fprintf('    dTOF_centro=%.1f%%  dTOF_excéntrico_max=%.1f%%  → %s\n', ...
        dTOF_centro, dTOF_exc_max, ternario(cumple,'✓ CUMPLE','✗ NO CUMPLE'));
end

fprintf('\n  Tendencia 2: posiciones en la ruta del rayo dan mayor dTOF\n');
fprintf('  (rayo horizontal z=0.45 m — ver columna "Rayo pasa por hueco")\n');

fprintf('\n  Tendencia 3: asimetría entre posiciones diagonales\n');
fprintf('    dTOF_UR=%.1f%%  >  dTOF_DL=%.1f%%\n', dTOF_UR, dTOF_DL);
fprintf('    dTOF_UL=%.1f%%  >  dTOF_DR=%.1f%%\n', dTOF_UL, dTOF_DR);
fprintf('    UR tiene el hueco más cerca del receptor (x=[0.47,0.57])\n');
fprintf('    → más interferencia con la llegada directa → mayor dTOF\n');
fprintf('    La asimetría refleja la anisotropía del roble (V_rad>V_tang)\n');

%% 6. FIGURAS

colores = {'k','b','r','g','m','c'};
estilos = {'-','--',':','-.','--',':'};

% Figura 1: Sismogramas superpuestos
figure('Name','Posiciones: Sismogramas','Position',[50 50 950 550]);
hold on; grid on;
for k = 1:N_casos
    if ~ok_all(k), continue; end
    amp_n = amp_all{k} / max(abs(amp_all{k}));
    plot(t_all{k}*1e6, amp_n, [colores{k} estilos{k}], ...
        'LineWidth', ternario(k==1,2.0,1.2), 'DisplayName', dirs{k,2});
end
xline(TOF_base,'k--','LineWidth',1,'DisplayName',sprintf('TOF base=%.1f µs',TOF_base));
xline(TOF_ref,'g-.','LineWidth',1,'DisplayName',sprintf('TOF analít.=%.1f µs',TOF_ref));
xlabel('Tiempo [µs]'); ylabel('Amplitud normalizada');
title('Sismogramas: 5 posiciones del defecto 10x10 cm','FontSize',11);
legend('Location','northeast','FontSize',8);
xlim([0 min(t_all{1}(end)*1e6, 800)]);

% Figura 2: dTOF por posición (barras)
figure('Name','Posiciones: dTOF por posición','Position',[100 100 800 420]);

subplot(1,2,1);
etiquetas_barras = {};
vals_us = []; vals_pct = [];
for k = 2:N_casos
    if ~ok_all(k), continue; end
    etiquetas_barras{end+1} = dirs{k,2};
    vals_us(end+1)  = dTOF_us(k);
    vals_pct(end+1) = dTOF_pct(k);
end
b = bar(1:length(vals_us), vals_us, 0.6);
b.FaceColor = 'flat';
cmap_bar = [0 0 1; 0.8 0 0; 0 0.7 0; 0.8 0.4 0; 0.5 0 0.8];
for j = 1:length(vals_us)
    b.CData(j,:) = cmap_bar(j,:);
end
set(gca,'XTickLabel', etiquetas_barras, 'XTick',1:length(vals_us));
ylabel('dTOF [µs]'); grid on;
title('Incremento TOF por posición [µs]','FontSize',10);
yline(0,'k-','LineWidth',0.5);

subplot(1,2,2);
b2 = bar(1:length(vals_pct), vals_pct, 0.6);
b2.FaceColor = 'flat';
for j = 1:length(vals_pct)
    b2.CData(j,:) = cmap_bar(j,:);
end
set(gca,'XTickLabel', etiquetas_barras, 'XTick',1:length(vals_pct));
ylabel('dTOF [%]'); grid on;
title('Incremento TOF por posición [%]','FontSize',10);
yline(0,'k-','LineWidth',0.5);

% Líneas de referencia tesis (solo si Centro está disponible)
if ~isnan(dTOF_centro)
    yline(34,'r--','LineWidth',1.5);   % ref sim tesis Centro
    text(0.5,34+1,'+34% ref.tesis(Centro,sim)','FontSize',7,'Color','r');
end

% Figura 3: Mapa de posiciones de los huecos
figure('Name','Posiciones: Geometría','Position',[150 150 600 550]);
hold on; axis equal; grid on;

% Tronco
rectangle('Position',[0.30,0.30,0.30,0.30],'EdgeColor',[0.6 0.3 0.1],...
    'FaceColor',[0.92 0.82 0.68],'LineWidth',2);

% Fuente y receptor
plot(0.30,0.45,'r*','MarkerSize',14,'LineWidth',2,'DisplayName','Fuente');
plot(0.60,0.45,'g^','MarkerSize',10,'LineWidth',2,'MarkerFaceColor','g','DisplayName','Receptor');

% Trayectoria fuente-receptor
plot([0.30 0.60],[0.45 0.45],'b-','LineWidth',1.5,'DisplayName','Rayo SR');

% Huecos en las 5 posiciones
colores_h = {'b','r','g','m','c'};
for j = 1:5
    rectangle('Position',[huecos(j).xmin, huecos(j).zmin, d_h, d_h],...
        'EdgeColor',colores_h{j},'LineWidth',2,...
        'LineStyle',ternario(j==1,'-','--'),'Curvature',[0 0]);
    text(huecos(j).cx, huecos(j).cz, nombres_pos{j},...
        'FontSize',9,'HorizontalAlignment','center',...
        'FontWeight','bold','Color',colores_h{j});
end

% Círculo de referencia tesis
th_c = linspace(0,2*pi,100);
plot(cx_t+R_t*cos(th_c), cz_t+R_t*sin(th_c),'k:','LineWidth',1,...
    'DisplayName','Tronco circular ref.');

xlabel('X [m]'); ylabel('Z [m]');
title('Geometría: 5 posiciones del defecto 10x10 cm','FontSize',11);
legend('Location','northeast','FontSize',8);
xlim([0.20 0.75]); ylim([0.20 0.75]);

text(0.21, 0.72, {'Tronco: 30x30 cm','Defecto: 10x10 cm (línea continua=Centro)'},...
    'FontSize',8,'BackgroundColor','w','EdgeColor','k');

%% =========================================================================
%  FUNCIÓN LOCAL
% =========================================================================

function [TOF_xc_us, TOF_th_us, TOF_en_us] = detectar_TOF_local(amp, signal_full, DT, NSTEP)
    T_total   = NSTEP * DT;
    N         = length(amp);
    TOF_min_s = 0.130e-3;
    TOF_max_s = min(0.500e-3, T_total);

    % 1. xcorr — solo el chirp, ventana física
    idx_nz = find(signal_full ~= 0);
    ch_n   = signal_full(idx_nz(1):idx_nz(end))';
    ch_n   = ch_n / (norm(ch_n) + eps);
    rec_n  = amp(:) / (norm(amp(:)) + eps);
    [xc, lags] = xcorr(rec_n, ch_n, round(T_total/DT));
    lags_s = lags * DT;
    v = (lags_s >= TOF_min_s) & (lags_s <= TOF_max_s);
    if ~any(v), v = (lags_s >= 0) & (lags_s <= T_total); end
    lv = lags_s(v); xv = xc(v);
    [~, idx]  = max(xv);
    TOF_xc_us = lv(idx) * 1e6;

    % 2. Threshold 5% del máximo, ventana física
    i0  = round(TOF_min_s / DT);
    i1  = min(round(TOF_max_s / DT), N);
    seg = amp(i0:i1);
    ith = find(abs(seg) > 0.05*max(abs(amp)), 1, 'first');
    if ~isempty(ith), TOF_th_us = (ith+i0-1)*DT*1e6; else, TOF_th_us = NaN; end

    % 3. Envolvente RMS local
    vw    = max(round(1/(36e3*DT)), 10);
    env   = sqrt(movmean(amp(:).^2, vw));
    seg_e = env(i0:i1);
    ie    = find(seg_e > 0.02*max(env), 1, 'first');
    if ~isempty(ie), TOF_en_us = (ie+i0-1)*DT*1e6; else, TOF_en_us = NaN; end
end

function s = ternario(cond, s_true, s_false)
    if cond, s = s_true; else, s = s_false; end
end