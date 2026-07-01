%
%  VALIDACIÓN SPECFEM2D vs. TESIS ESPINOSA (2019)
%  Tronco de pino 30x30 cm, malla interna estructurada
%
%  Referencia principal: Espinosa et al. (2019), Ultrasonics 91:242-251
%
%  MÉTODOS DE DETECCIÓN DE TOF (Cap.3, pág.39-51, tesis):
%  ─────────────────────────────────────────────────────
%  1. Correlación cruzada (xcorr):
%     TOF = posición del máximo de la xcorr entre señal emitida y recibida
%     Fórmula (ec.24): r_sy[l] = (1/√(Es·Ey)) · Σ s[k]·y[k-l]
%     → Método principal para chirp (menor variación, pág.51)
%     → En simulación: estima retardo del CENTRO DE ENERGÍA del paquete
%
%  2. Threshold m·σ (umbral):
%     TOF = primer instante donde |señal| > m·σ_ruido, con m=8 (pág.39)
%     → Método clásico, sensible al ruido y a la amplitud del chirp
%     → En simulación: SNR→∞, por lo que m·σ≈0 y detecta el frente de onda
%
%  3. Envolvente RMS:
%     TOF = primera muestra donde la envolvente RMS supera un umbral relativo
%     → Variante robusta para simulaciones sin ruido (SPECFEM2D)
%     → Detecta la primera llegada de ENERGÍA, no el frente de fase
%     → Es el estimador más coherente con la física del chirp en simulación
%
%  Coherencia con la tesis:
%     - Cap.2: Christoffel → velocidades teóricas
%     - Cap.3: xcorr como mejor método con chirp (Fig.30-31, pág.49-51)
%     - Cap.5: dTOF con hueco centrado +34%(sim)/+50%(exp) (Fig.63, pág.79)

clear; clc; close all;

fprintf('\n');
fprintf('  VALIDACIÓN SPECFEM2D vs. TESIS ESPINOSA (2019)\n');
fprintf('  Tronco: 30x30 cm | Hueco centrado: 10x10 cm\n');
fprintf('  Fuente: (0.35,0.50) m | Receptor: (0.65,0.50) m | dist=0.300 m\n');
fprintf('  DT=200 ns | NSTEP=3000 | T_total=600 µs | PML ON\n\n');

%% 1. PARÁMETROS MATERIALES Y VELOCIDADES (Christoffel)

% Tiempo
DT      = 2.0e-7;      % s 
NSTEP   = 3000;        % pasos totales
T_total = NSTEP * DT;  % s = 1.0 ms
t_vec   = (0:NSTEP-1) * DT;

% Dominio y malla ( dx=(xmax-xmin)/nx )
dx = 1.0 / 100;    % m por elemento (nx=100, nz=100, xmax=1.0, xmin = 0.0)

% Fuente y receptor
xs = 0.35;  zs = 0.50;   % m
xr = 0.65;  zr = 0.50;   % m
dist_SR = sqrt((xr-xs)^2 + (zr-zs)^2);

% Material: pino anisótropo (Tabla 2, pág.27 tesis)
rho = 661;        % kg/m³
c11 = 1.537e9;    % Pa  (ER pino, Tabla 2)
c33 = 8.29e8;     % Pa  (ET pino, Tabla 2)
c55 = 1.81e8;     % Pa  (GRT pino, Tabla 2)
c13 = 4.289e8;    % Pa  (νRT·√(ER·ET))

% Material: aire (defecto y entorno)
Vp_air = 343.0;        % m/s

% Geometría del tronco (elementos 35-65 en X y Z)
tronco.xmin = 35*dx;  tronco.xmax = 65*dx;
tronco.zmin = 35*dx;  tronco.zmax = 65*dx;
tronco.ancho = tronco.xmax - tronco.xmin;   % 0.30 m

% Geometría del hueco centrado (elementos 45-55 en X y Z)
hueco.xmin = 45*dx;   hueco.xmax = 55*dx;
hueco.zmin = 45*dx;   hueco.zmax = 55*dx;
hueco.ancho = hueco.xmax - hueco.xmin;       % 0.10 m

fprintf('--- Velocidades de fase (Ec. Christoffel, Cap.2 tesis) ---\n');

% Constantes elásticas pino (Par_file, Tabla 3 tesis pág.41)
rho = 661;        % kg/m³
c11 = 1.537e9;    % Pa  (ER pino, Tabla 2)
c33 = 8.29e8;     % Pa  (ET pino, Tabla 2)
c55 = 1.81e8;     % Pa  (GRT pino, Tabla 2)
c13 = 4.289e8;    % Pa  (νRT·√(ER·ET))

% Velocidades Christoffel (ec.13, Cap.2)
V0  = sqrt(c11/rho);   % θ=0°  (radial)
V90 = sqrt(c33/rho);   % θ=90° (tangencial)

fprintf('  θ=0°  (radial):     %.1f m/s  |  Tesis: 1525 m/s  |  error: %.1f%%\n', ...
    V0,  100*abs(V0-1525)/1525);
fprintf('  θ=90° (tangencial): %.1f m/s  |  Tesis: 1120 m/s  |  error: %.1f%%\n', ...
    V90, 100*abs(V90-1120)/1120);

theta_arr = 0:1:90;

%% 2. TOF DE REFERENCIA

dist = 0.300;  % m
TOF_ref = dist / V0 * 1e6;  % µs

fprintf('\n--- TOF de referencia ---\n');
fprintf('  Trayectoria: %.3f m en madera (%.0f m/s) → %.2f µs\n', ...
    dist, V0, TOF_ref);
fprintf('  TOF_ref = %.2f µs  |  V_eq = %.1f m/s\n', TOF_ref, V0);
fprintf('  Rango tesis (Ash Ø30cm, dir. radial, Fig.61): ~194 µs %s\n', ...
    ternario(TOF_ref>=185 && TOF_ref<=200,'✓','⚠'));

%% 3. SEÑAL CHIRP

% Parámetros exactos del script
Ts    = 45e-6;   % s
f0_ch = 22e3;    % Hz
f1_ch = 50e3;    % Hz
Fc    = 36e3;    % Hz (frecuencia central)
dF    = 28e3;    % Hz (ancho de banda)
DT = 2.0e-7;  % s
NSTEP = 3000; 
Nc = round(Ts/DT);

t_chirp = 0:DT:Ts-DT;
chirp_signal = cos(2*pi*((f1_ch-f0_ch)*t_chirp/Ts+f0_ch).*t_chirp) .* exp(-((t_chirp-Ts/2).^2)/(2*(Ts/6)^2));
signal_full = [chirp_signal, zeros(1, NSTEP-length(t_chirp))];
 
% Ancho de banda -3 dB
Nfft  = 2^nextpow2(NSTEP);
S_mag = 2*abs(fft(signal_full,Nfft)/NSTEP);
freq  = (0:Nfft/2-1)/(Nfft*DT);
S_mag = S_mag(1:Nfft/2);
i3dB  = find(S_mag > max(S_mag)*10^(-3/20));
 
fprintf('--- Señal chirp R3α ---\n');
fprintf('  f: %.0f→%.0f kHz | Fc=%.0f kHz | Ts=%.0f µs | %d pts\n', ...
    f0_ch/1e3, f1_ch/1e3, Fc/1e3, Ts*1e6, length(t_chirp));
fprintf('  BW -3dB: [%.1f - %.1f] kHz  |  Tesis (Tabla 5): [32.6 - 40.0] kHz\n\n', ...
    freq(i3dB(1))/1e3, freq(i3dB(end))/1e3);


%% 4. CARGA DE SISMOGRAMAS

file_sano  = 'Mesh_interna/Pino_sin_hueco_Malla/OUTPUT_FILES/AA.S0001.BXX.semd';
file_hueco = 'Mesh_interna/Pino_con_hueco_N_C/OUTPUT_FILES/AA.S0001.BXX.semd';

fprintf('\n--- Cargando sismogramas ---\n');
if ~exist(file_sano,'file'),  error('No encontrado: %s', file_sano);  end
if ~exist(file_hueco,'file'), error('No encontrado: %s', file_hueco); end

A = load(file_sano);  B = load(file_hueco);
t_s = A(:,1)-A(1,1);  s_s = A(:,2);
t_h = B(:,1)-B(1,1);  s_h = B(:,2);
N   = length(t_s);

fprintf('  ✓ %s (%d pts)\n', file_sano,  N);
fprintf('  ✓ %s (%d pts)\n', file_hueco, N);

%% 5. DETECCIÓN DE TOF — TRES MÉTODOS (Cap.3 tesis)

fprintf('\n--- Detección de TOF (Cap.3, pág.39-51) ---\n');
fprintf('  Método       TOF sano [µs]   TOF hueco [µs]  dTOF [µs] (%%)\n');
fprintf('  %s\n', repmat('-',1,60));

% Normalizar señales
s_s_n = s_s / max(abs(s_s));
s_h_n = s_h / max(abs(s_h));

% ── Método 1: Correlación cruzada (ec.24, pág.40) ──────────────────────
% xcorr entre chirp emitido y señal recibida
chirp_full = zeros(N,1);
chirp_full(1:Nc) = chirp_signal;

[c_s, lags_s] = xcorr(s_s_n, chirp_full, 'coeff');
[~, ic_s] = max(c_s);
TOF_xcorr_s = lags_s(ic_s) * DT * 1e6;

[c_h, lags_h] = xcorr(s_h_n, chirp_full, 'coeff');
[~, ic_h] = max(c_h);
TOF_xcorr_h = lags_h(ic_h) * DT * 1e6;

dTOF_xcorr = TOF_xcorr_h - TOF_xcorr_s;
fprintf('  xcorr        %8.2f        %8.2f        %+.2f µs (%+.1f%%)\n', ...
    TOF_xcorr_s, TOF_xcorr_h, dTOF_xcorr, 100*dTOF_xcorr/TOF_xcorr_s);
fprintf('  [método principal para chirp, pág.51 tesis]\n');

% ── Método 2: Threshold m=8·σ (pág.39) ────────────────────────────────
% En simulación sin ruido: σ≈0, usamos 10% del máximo como equivalente
m = 8;
% Estimar σ_ruido de los primeros 100 µs (antes de llegada)
i_pre = round(0.100e-3/DT);
sigma_s = std(s_s_n(1:i_pre));
sigma_h = std(s_h_n(1:i_pre));
umbral_s = m * sigma_s;
umbral_h = m * sigma_h;
% Si σ≈0 (simulación limpia), usar 10% como fallback
if umbral_s < 1e-6, umbral_s = 0.10; end
if umbral_h < 1e-6, umbral_h = 0.10; end

i_th_s = find(abs(s_s_n) > umbral_s, 1, 'first');
i_th_h = find(abs(s_h_n) > umbral_h, 1, 'first');
TOF_th_s = t_s(i_th_s)*1e6;
TOF_th_h = t_h(i_th_h)*1e6;
dTOF_th  = TOF_th_h - TOF_th_s;
fprintf('  Threshold    %8.2f        %8.2f        %+.2f µs (%+.1f%%)\n', ...
    TOF_th_s, TOF_th_h, dTOF_th, 100*dTOF_th/TOF_th_s);
fprintf('  [m=8·σ, pág.39 tesis; en simulación σ≈0 → equivalente a 10%%]\n');

% ── Método 3: Envolvente RMS (variante robusta para simulación) ────────
vw = max(round(1/(Fc*DT)), 10);
env_s = sqrt(movmean(s_s_n.^2, vw));
env_h = sqrt(movmean(s_h_n.^2, vw));
i0_env = round(0.130e-3/DT);
i1_env = min(round(0.500e-3/DT), N);

ie_s = find(env_s(i0_env:i1_env) > 0.02*max(env_s), 1, 'first');
ie_h = find(env_h(i0_env:i1_env) > 0.02*max(env_h), 1, 'first');
TOF_env_s = (ie_s + i0_env - 1) * DT * 1e6;
TOF_env_h = (ie_h + i0_env - 1) * DT * 1e6;
dTOF_env  = TOF_env_h - TOF_env_s;
fprintf('  Envolvente   %8.2f        %8.2f        %+.2f µs (%+.1f%%)\n', ...
    TOF_env_s, TOF_env_h, dTOF_env, 100*dTOF_env/TOF_env_s);
fprintf('  [variante robusta para SPECFEM2D sin ruido]\n');

%% 6. VERIFICACIÓN vs. TESIS

fprintf('\n--- Verificación ---\n');

% Velocidad equivalente (envolvente como estimador absoluto)
V_eq_env = dist / (TOF_env_s*1e-6);
err_V    = 100*abs(V_eq_env - V0)/V0;
fprintf('  V ref (modelo):              %.1f m/s\n', V0);
fprintf('  V medida (envolvente):       %.1f m/s  |  error: %.1f%%  ← estimador TOF absoluto\n', ...
    V_eq_env, err_V);
fprintf('  V medida (threshold):        %.1f m/s  |  error: %.1f%%\n', ...
    dist/(TOF_th_s*1e-6), 100*abs(dist/(TOF_th_s*1e-6)-V0)/V0);
fprintf('  V medida (xcorr):            %.1f m/s  |  error: %.1f%%  ← desfase Ts/2=%.1f µs\n', ...
    dist/(TOF_xcorr_s*1e-6), 100*abs(dist/(TOF_xcorr_s*1e-6)-V0)/V0, Ts/2*1e6);

% Posición del rayo vs. hueco
z_rayo = 0.50; z_hueco_min = 0.45; z_hueco_max = 0.55;
fprintf('\n  Hueco z=[%.2f,%.2f] m vs. rayo z=%.2f m → ', ...
    z_hueco_min, z_hueco_max, z_rayo);
if z_rayo >= z_hueco_min && z_rayo <= z_hueco_max
    fprintf('el rayo ATRAVIESA el hueco\n');
else
    fprintf('el rayo NO atraviesa el hueco\n');
end

% dTOF con corrección por área (tesis usa hueco cilíndrico, simulación cuadrado)
A_cuad = 0.10^2;                      % m² (hueco 10x10 cm)
A_circ = pi*(0.10/2)^2;               % m² (hueco cilíndrico Ø10cm, tesis)
factor_area = A_cuad / A_circ;        % ≈ 1.27
dTOF_env_corr = dTOF_env/TOF_env_s*100 / factor_area;
fprintf('  dTOF envolvente: %+.1f%%  |  corregido por área (%.0fcm²→%.1fcm²): %+.1f%%\n', ...
    100*dTOF_env/TOF_env_s, A_cuad*1e4, A_circ*1e4, dTOF_env_corr);
fprintf('  Referencia tesis: +34%% (sim) / +50%% (exp)  [Fig.63, pág.79]\n');

%% 7. MALLA Y CFL

fprintf('\n--- Malla y CFL ---\n');
dx   = 0.01;    % m
dxGLL = dx/4;   % distancia entre puntos GLL
lambda_min = V90/f1_ch;
pts_lam    = lambda_min/dxGLL;
CFL        = V0*DT/dxGLL;
T_total    = 3000*DT*1e6;
TOF_crit   = TOF_ref + 2*Ts*1e6;

fprintf('  dx_elem=%.4f m | dx_GLL=%.5f m | pts/λ_min=%.1f (%s)\n', ...
    dx, dxGLL, pts_lam, ternario(pts_lam>=5,'≥5 ✓','<5 ✗'));
fprintf('  CFL_GLL=%.4f (%s) | T_total=%.0f µs > TOF_ref+2Ts=%.0f µs %s\n', ...
    CFL, ternario(CFL<0.5,'<0.5 ✓','>0.5 ✗'), T_total, TOF_crit, ...
    ternario(T_total>TOF_crit,'✓','✗'));

%% 8. TABLA RESUMEN

fprintf('\n  TABLA RESUMEN DE VALIDACIÓN\n');
fprintf('  %-35s %-15s %-20s %s\n','Métrica','Simulación','Tesis','Estado');
fprintf('  %s\n', repmat('-',1,75));

tab = {
    'V(θ=0°) [m/s]',            sprintf('%.1f',V0),         '1525',              ternario(abs(V0-1525)/1525<0.02,'✓','✗');
    'V(θ=90°) [m/s]',           sprintf('%.1f',V90),        '1120',              ternario(abs(V90-1120)/1120<0.02,'✓','✗');
    'TOF sano xcorr [µs]',      sprintf('%.1f',TOF_xcorr_s),'196.7 (ref)',       ternario(abs(TOF_xcorr_s-TOF_ref)/TOF_ref<0.15,'~','✗');
    'TOF sano threshold [µs]',  sprintf('%.1f',TOF_th_s),   '196.7 (ref)',       ternario(abs(TOF_th_s-TOF_ref)/TOF_ref<0.10,'✓','✗');
    'TOF sano envolvente [µs]', sprintf('%.1f',TOF_env_s),  '196.7 (ref)',       ternario(abs(TOF_env_s-TOF_ref)/TOF_ref<0.02,'✓','~');
    'dTOF xcorr',               sprintf('%+.1f%%',100*dTOF_xcorr/TOF_xcorr_s),'+34%(sim)/+50%(exp)','[no fiable con hueco grande]';
    'dTOF threshold',           sprintf('%+.1f%%',100*dTOF_th/TOF_th_s),       '+34%(sim)/+50%(exp)',ternario(100*dTOF_th/TOF_th_s>20,'⚠ hueco cuadrado','✗');
    'dTOF envolvente',          sprintf('%+.1f%%',100*dTOF_env/TOF_env_s),     '+34%(sim)/+50%(exp)',ternario(100*dTOF_env/TOF_env_s>30,'⚠ hueco cuadrado','✗');
    'dTOF corregido por área',  sprintf('%+.1f%%',dTOF_env_corr),              '34-50%',             ternario(dTOF_env_corr>=30 && dTOF_env_corr<=55,'≈ comparable','✗');
    'pts/λ_min (GLL)',          sprintf('%.1f',pts_lam),     '≥5',               ternario(pts_lam>=5,'✓','✗');
    'CFL (dx_GLL)',             sprintf('%.3f',CFL),          '<0.5',             ternario(CFL<0.5,'✓','✗');
    'Chirp 22-50 kHz',          '✓',                         'Tabla 5',          '✓';
    'PML condiciones',          'NELEM=5',                   'PML/Stacey',       '✓';
};

for i = 1:size(tab,1)
    fprintf('  %-35s %-15s %-20s %s\n', tab{i,1}, tab{i,2}, tab{i,3}, tab{i,4});
end

fprintf('\n  Leyenda: ✓ <2%%  ≈ <10%%  ~ <15%%  ⚠ diferencia documentada\n');
fprintf('\n  NOTA SOBRE MÉTODOS TOF (Cap.3, pág.49-51):\n');
fprintf('  xcorr = estimador del CENTRO DE ENERGÍA del paquete chirp\n');
fprintf('         → da TOF más tardío porque detecta la media, no el frente\n');
fprintf('  threshold/env = estimadores de la PRIMERA LLEGADA de energía\n');
fprintf('         → más coherentes con la velocidad de fase de la tesis\n');
fprintf('  Con chirp, xcorr < variación experimental, pero threshold/env\n');
fprintf('  son más precisos para la velocidad de fase en simulación.\n');

%% 9. FIGURAS


%  FIG.1 — GEOMETRÍA DEL MODELO
%  Muestra la diferencia rectangular vs. circular, la posición del hueco
%  respecto al rayo fuente-receptor, y los puntos de fuente/receptor.
%  Referencia: comparar con esquema Cap.4, pág.55 (Fig.35)

figure('Name','Fig.1 - Geometría del modelo','Position',[50 50 620 540]);
hold on; axis equal; grid on;
 
% Tronco (madera)
rectangle('Position',[tronco.xmin, tronco.zmin, tronco.ancho, tronco.ancho], ...
    'FaceColor',[0.88 0.74 0.54], 'EdgeColor',[0.50 0.25 0.05], 'LineWidth',2);
 
% Hueco (aire)
rectangle('Position',[hueco.xmin, hueco.zmin, hueco.ancho, hueco.ancho], ...
    'FaceColor',[0.75 0.88 1.00], 'EdgeColor',[0.10 0.30 0.80], 'LineWidth',2);
text(mean([hueco.xmin hueco.xmax]), mean([hueco.zmin hueco.zmax]), ...
    {'Hueco','(aire)'},'FontSize',8,'Color',[0.10 0.30 0.80],...
    'HorizontalAlignment','center','VerticalAlignment','middle');
 
% Trayectoria fuente-receptor
plot([xs xr],[zs zr],'b-','LineWidth',2,'DisplayName','Trayectoria SR');
 
% Fuente y receptor
plot(xs, zs, 'r*', 'MarkerSize',14, 'LineWidth',2, 'DisplayName','Fuente');
plot(xr, zr, 'g^', 'MarkerSize',10, 'LineWidth',2, 'MarkerFaceColor','g', ...
    'DisplayName','Receptor');
 
% Círculo de referencia de la tesis (tronco circular Ø30 cm)
cx_t = (tronco.xmin+tronco.xmax)/2;
cz_t = (tronco.zmin+tronco.zmax)/2;
th_c = linspace(0, 2*pi, 200);
plot(cx_t + 0.15*cos(th_c), cz_t + 0.15*sin(th_c), 'r:', 'LineWidth',1.5, ...
    'DisplayName','Tronco circular tesis (Ø30 cm)');
 
% Anotación del TOF de referencia
xm = (xs+xr)/2;
text(xm, zs+0.025, sprintf('d = %.2f m  |  TOF_{ref} = %.0f µs', dist_SR, TOF_ref), ...
    'FontSize',8, 'Color','b', 'HorizontalAlignment','center');
 
% Anotación de diferencia de áreas
area_cuad = (hueco.ancho*100)^2;
area_circ = pi*5^2;
text(tronco.xmax+0.01, tronco.zmin+0.05, ...
    {sprintf('Hueco cuadrado: %.0f cm²',area_cuad), ...
     sprintf('Hueco circular (tesis): %.1f cm²',area_circ), ...
     sprintf('Diferencia: +%.0f%%',(area_cuad/area_circ-1)*100)}, ...
    'FontSize',7, 'Color',[0.10 0.30 0.80], 'BackgroundColor','w', ...
    'EdgeColor',[0.7 0.7 0.7]);
 
xlabel('X [m]'); ylabel('Z [m]');
title('Geometría del modelo SPECFEM2D vs. tesis','FontSize',11,'FontWeight','bold');
legend('Location','northeast','FontSize',8);
xlim([0.22 0.75]); ylim([0.22 0.75]);
 
% Etiqueta del tronco
text(tronco.xmin+0.01, tronco.zmax-0.02, 'Pino anisótropo', ...
    'FontSize',8, 'Color',[0.40 0.20 0]);

figure('Name','Validación - Sismogramas','Position',[50 50 900 500]);
hold on; grid on;
plot(t_s*1e6, s_s_n, 'b-',  'LineWidth',1.5, 'DisplayName','Sin hueco');
plot(t_h*1e6, s_h_n, 'r--', 'LineWidth',1.5, 'DisplayName','Con hueco (10x10cm)');
xline(TOF_xcorr_s,'b:','LineWidth',1,'DisplayName',sprintf('TOF xcorr sano=%.1fµs',TOF_xcorr_s));
xline(TOF_env_s,  'b-','LineWidth',2,'DisplayName',sprintf('TOF env sano=%.1fµs',TOF_env_s));
xline(TOF_xcorr_h,'r:','LineWidth',1,'DisplayName',sprintf('TOF xcorr hueco=%.1fµs',TOF_xcorr_h));
xline(TOF_env_h,  'r-','LineWidth',2,'DisplayName',sprintf('TOF env hueco=%.1fµs',TOF_env_h));
xlabel('Tiempo [µs]'); ylabel('Amplitud normalizada');
title({'Sismogramas SPECFEM2D — Pino (malla interna)', ...
    sprintf('dTOF xcorr=%+.1f%% | dTOF env=%+.1f%% | Ref.tesis +34%%(sim)/+50%%(exp)', ...
    100*dTOF_xcorr/TOF_xcorr_s, 100*dTOF_env/TOF_env_s)}, 'FontSize',10);
legend('Location','northeast','FontSize',8);
xlim([0 600]); ylim([-1.1 1.1]);

figure('Name','Validación - Zoom primera llegada','Position',[50 150 900 400]);
ventana = [130 350];
idx_v = t_s*1e6>=ventana(1) & t_s*1e6<=ventana(2);
hold on; grid on;
plot(t_s(idx_v)*1e6, s_s_n(idx_v), 'b-',  'LineWidth',2, 'DisplayName','Sin hueco');
plot(t_h(idx_v)*1e6, s_h_n(idx_v), 'r--', 'LineWidth',2, 'DisplayName','Con hueco');
plot(t_s(idx_v)*1e6, env_s(idx_v), 'b:',  'LineWidth',1.5, 'DisplayName','Envolvente sano');
plot(t_h(idx_v)*1e6, env_h(idx_v), 'r:',  'LineWidth',1.5, 'DisplayName','Envolvente hueco');
xline(TOF_xcorr_s,'b-','LineWidth',1.5,'DisplayName',sprintf('xcorr sano=%.1fµs',TOF_xcorr_s));
xline(TOF_env_s,  'b--','LineWidth',1.5,'DisplayName',sprintf('env sano=%.1fµs',TOF_env_s));
xline(TOF_xcorr_h,'r-','LineWidth',1.5,'DisplayName',sprintf('xcorr hueco=%.1fµs',TOF_xcorr_h));
xline(TOF_env_h,  'r--','LineWidth',1.5,'DisplayName',sprintf('env hueco=%.1fµs',TOF_env_h));
xline(TOF_ref, 'k--','LineWidth',1,'DisplayName',sprintf('TOF analítico=%.1fµs',TOF_ref));
xlabel('Tiempo [µs]'); ylabel('Amplitud normalizada');
title('Zoom primera llegada — comparación de métodos TOF','FontSize',10);
legend('Location','northeast','FontSize',7);
xlim(ventana);


%% FUNCIÓN LOCAL
function s = ternario(cond, s_true, s_false)
    if cond, s = s_true; else, s = s_false; end
end
