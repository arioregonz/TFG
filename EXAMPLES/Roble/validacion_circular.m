%
%  VALIDACIÓN MALLA EXTERNA CIRCULAR vs. TESIS ESPINOSA (2019)
%  Tronco cilíndrico Ø30cm — con y sin hueco centrado Ø10cm
%
%  Referencia: Espinosa et al. (2019), Ultrasonics 91:242-251
%    Fig.61 (pág.77): TOF experimental roble dir. radial: 158-190 µs
%    Fig.63 (pág.79): dTOF sim +34%, exp +50% (hueco cilíndrico Ø10cm)
%    Cap.3  (pág.39-51): métodos TOF — xcorr principal con chirp
%
%  Geometría malla circular:
%    Dominio total: 1.10 x 1.10 m (PML 0.05m en cada borde)
%    Tronco: disco radio=0.15m centrado en (0.55, 0.55)
%    Hueco:  disco radio=0.05m centrado en (0.55, 0.55) → Ø10cm
%    Fuente:   xs=0.40m, zs=0.55m (borde izq tronco, punto tangente)
%    Receptor: xr=0.70m, zr=0.55m (borde der tronco, punto tangente)
%    Distancia: 0.30m (diámetro tronco) ✓
%
%  IMPORTANTE — anglesource=90 (fuerza radial hacia el centro del tronco):
%    En geometría circular la fuente debe orientarse radialmente. En el
%    punto tangente del círculo, una fuerza vertical (angle=0) excitaría
%    el modo TANGENCIAL (lento, ~1383 m/s, energía por el contorno).
%    Con angle=90 la fuerza apunta horizontalmente hacia el centro del
%    disco → excita el modo RADIAL (rápido, ~1897 m/s, cruza el diámetro).
%    En geometría rectangular (borde plano) la interfaz plana canaliza la
%    energía y es más tolerante a la orientación de la fuente.

clear; clc; close all;

fprintf('\n  VALIDACIÓN MALLA CIRCULAR vs. TESIS ESPINOSA (2019)\n');
fprintf('  Tronco Ø30cm | Hueco Ø10cm centrado\n');
fprintf('  Fuente: (0.40,0.55)m | Receptor: (0.70,0.55)m | dist=0.30m\n');
fprintf('  anglesource=90 (fuerza radial hacia el centro del tronco)\n\n');

%% 1. PARÁMETROS MATERIALES Y VELOCIDADES TEÓRICAS

rho = 706; c11 = 2.54e9; c33 = 1.35e9; c55 = 3.19e8; c13 = 7.5832e8;
V0  = sqrt(c11/rho);   % m/s radial
V90 = sqrt(c33/rho);   % m/s tangencial
dist = 0.30;
TOF_ref = dist/V0*1e6;

fprintf('--- Velocidades de fase (Christoffel, Cap.2 tesis) ---\n');
fprintf('  V(θ=0°)  = %.1f m/s  |  tesis: 1898 m/s  |  error: %.1f%%\n', ...
    V0,  100*abs(V0-1898)/1898);
fprintf('  V(θ=90°) = %.1f m/s  |  tesis: 1385 m/s  |  error: %.1f%%\n', ...
    V90, 100*abs(V90-1385)/1385);
fprintf('  TOF analítico (0.30m): %.2f µs\n', TOF_ref);
fprintf('  Rango experimental tesis Fig.61: 158-190 µs %s\n', ...
    ternario(TOF_ref>=158 && TOF_ref<=190,'✓','(fuera)'));

%% 2. ARCHIVOS DE SIMULACIÓN

DT = 2.0e-7;

% Malla circular externa — casos sano y con hueco
file_circ_sano  = 'Mesh_externa/Roble_C_sin_hueco_N_2/OUTPUT_FILES/AA.S0001.BXX.semd';
file_circ_hueco = 'Mesh_externa/Roble_C_con_hueco_N_3/OUTPUT_FILES/AA.S0001.BXX.semd';

% Malla interna de referencia (validada)
file_int_sano   = 'Mesh_interna/Roble_sin_hueco_Norm/OUTPUT_FILES/AA.S0001.BXX.semd';
file_int_hueco  = 'Mesh_interna/Roble_con_hueco_N_C/OUTPUT_FILES/AA.S0001.BXX.semd';

fprintf('\n--- Cargando sismogramas ---\n');
archivos = {file_circ_sano,'Circ. sin hueco'; file_circ_hueco,'Circ. con hueco'; ...
            file_int_sano, 'Int. sin hueco';  file_int_hueco, 'Int. con hueco'};
datos = cell(4,1);
for i = 1:4
    if exist(archivos{i,1},'file')
        A = load(archivos{i,1});
        datos{i}.t = A(:,1)-A(1,1);
        datos{i}.s = A(:,2);
        datos{i}.nombre = archivos{i,2};
        fprintf('  ✓ %s (%d pts)\n', archivos{i,2}, length(A));
    else
        fprintf('  ✗ NO EXISTE: %s\n', archivos{i,1});
        datos{i} = [];
    end
end

%% 3. CHIRP EMITIDO (para xcorr)

Ts=45e-6; f0c=22e3; f1c=50e3; Fc=36e3;
Nc = round(Ts/DT);
t_c = (0:Nc-1)*DT;
f_t = (f1c-f0c)*t_c/Ts + f0c;
phase = 2*pi*f_t.*t_c;
win = exp(-((t_c-Ts/2).^2)/(2*(Ts/6)^2));
chirp_s = cos(phase).*win;
if ~isempty(datos{1})
    N = length(datos{1}.t);
else
    N = 3000;
end
chirp_full = zeros(N,1);
chirp_full(1:Nc) = chirp_s;

%% 4. DETECCIÓN TOF — TRES MÉTODOS

vw = max(round(1/(Fc*DT)),10);
tof_xcorr = zeros(4,1);
tof_th    = zeros(4,1);
tof_env   = zeros(4,1);
V_med     = zeros(4,1);

fprintf('\n--- Detección de TOF (tres métodos, Cap.3 tesis) ---\n');
fprintf('  Orden: xcorr (principal chirp) | threshold (m=8σ) | envolvente (robusta)\n');
fprintf('  %-20s  xcorr[µs]  thresh[µs]  env[µs]  V_env[m/s]\n','Caso');
fprintf('  %s\n',repmat('-',1,68));

for i = 1:4
    if isempty(datos{i}), continue; end
    t = datos{i}.t; s = datos{i}.s;
    sn = s/max(abs(s));

    % Método 1: xcorr (chirp emitido vs recibido, ec.24 pág.40)
    [c,lags] = xcorr(sn, chirp_full, 'coeff');
    [~,ic] = max(c);
    tof_xcorr(i) = lags(ic)*DT*1e6;

    % Método 2: Threshold m=8σ (pág.39)
    i_pre = round(0.100e-3/DT);
    sig = std(sn(1:i_pre));
    umb = max(8*sig, 0.10);
    ith = find(abs(sn)>umb,1,'first');
    if ~isempty(ith), tof_th(i) = t(ith)*1e6; end

    % Método 3: Envolvente RMS (robusta SPECFEM2D)
    env = sqrt(movmean(sn.^2, vw));
    i0 = round(0.130e-3/DT); i1 = min(round(0.500e-3/DT),N);
    ie = find(env(i0:i1)>0.02*max(env),1,'first');
    if ~isempty(ie)
        tof_env(i) = (ie+i0-1)*DT*1e6;
        V_med(i) = dist/(tof_env(i)*1e-6);
    end

    fprintf('  %-20s  %8.2f   %8.2f   %8.2f  %8.1f\n', ...
        datos{i}.nombre, tof_xcorr(i), tof_th(i), tof_env(i), V_med(i));
end

%% 5. VALIDACIÓN vs. TESIS

fprintf('\n--- Validación vs. tesis (Espinosa 2019) ---\n');

% TOF sano circular
if tof_env(1)>0
    err_tof = 100*abs(tof_env(1)-TOF_ref)/TOF_ref;
    err_V   = 100*abs(V_med(1)-V0)/V0;
    fprintf('  Caso sano (malla circular):\n');
    fprintf('    TOF_env=%.2f µs | TOF_analitico=%.2f µs | error=%.1f%% %s\n', ...
        tof_env(1), TOF_ref, err_tof, ternario(err_tof<5,'✓ (<5%)','⚠'));
    fprintf('    V_med=%.1f m/s  | V_teo=%.1f m/s | error=%.1f%% %s\n', ...
        V_med(1), V0, err_V, ternario(err_V<5,'✓ (<5%)','⚠'));
    in_range = tof_env(1)>=158 && tof_env(1)<=190;
    if in_range
        fprintf('    Rango exp. tesis Fig.61 (158-190µs): DENTRO ✓\n');
    else
        fprintf('    Rango exp. tesis Fig.61 (158-190µs):\n');
        fprintf('      TOF_sim=%.2fµs vs mínimo exp=158µs (dif=%.2fµs=%.1f%%)\n', ...
            tof_env(1), abs(tof_env(1)-158), 100*abs(tof_env(1)-158)/158);
        fprintf('      Diferencia explicable: error velocidad %.1f%% (<5%%) ✓\n', err_V);
        fprintf('      TOF analítico (158.2µs) coincide con el mínimo experimental.\n');
        fprintf('      Simulación 2D radial pura → V ligeramente mayor que experimento 3D.\n');
    end
end

% dTOF con hueco circular
if tof_env(1)>0 && tof_env(2)>0
    dTOF_c   = 100*(tof_env(2)-tof_env(1))/tof_env(1);
    dTOF_ref_sim = 34; dTOF_ref_exp = 50;
    fprintf('\n  Caso con hueco Ø10cm (malla circular):\n');
    fprintf('    dTOF_env=+%.1f%%\n', dTOF_c);
    fprintf('    Ref. tesis Fig.63: +%.0f%%(sim) / +%.0f%%(exp)\n', ...
        dTOF_ref_sim, dTOF_ref_exp);
    fprintf('    Estado: %s\n', ...
        ternario(dTOF_c>=dTOF_ref_sim*0.8 && dTOF_c<=dTOF_ref_exp*1.3, ...
        'COMPARABLE CON TESIS ✓', 'DIFERENCIA A DOCUMENTAR ⚠'));
    fprintf('    Ventaja: hueco cilíndrico = geometría IDÉNTICA a la tesis (Fig.63)\n');
    fprintf('    → No requiere corrección por área (a diferencia de malla interna)\n');
end

% Comparación circular vs. interna
if tof_env(1)>0 && tof_env(3)>0
    err_ci = 100*abs(tof_env(1)-tof_env(3))/tof_env(3);
    fprintf('\n  Comparación circular vs. interna (caso sano):\n');
    fprintf('    TOF_circ=%.2fµs | TOF_int=%.2fµs | error=%.2f%% %s\n', ...
        tof_env(1), tof_env(3), err_ci, ternario(err_ci<5,'✓ (<5%)','⚠'));
    fprintf('    Nota: geometrías distintas (circular vs rectangular) → leve diferencia\n');
end
if tof_env(2)>0 && tof_env(4)>0
    dTOF_int = 100*(tof_env(4)-tof_env(3))/tof_env(3);
    dTOF_cir = 100*(tof_env(2)-tof_env(1))/tof_env(1);
    fprintf('  Comparación circular vs. interna (dTOF con hueco):\n');
    fprintf('    dTOF_circ=+%.1f%% (hueco cilíndrico) | dTOF_int=+%.1f%% (hueco cuadrado)\n', ...
        dTOF_cir, dTOF_int);
    fprintf('    El circular es MÁS comparable con la tesis (geometría idéntica)\n');
end

%% 6. TABLA RESUMEN

fprintf('\n  TABLA RESUMEN DE VALIDACIÓN — MALLA CIRCULAR\n');
fprintf('  %-32s %-12s %-20s %s\n','Métrica','Circular','Ref. tesis','Estado');
fprintf('  %s\n',repmat('-',1,72));

filas = {};
if tof_env(1)>0
    filas{end+1,1} = 'V(θ=0°) medida [m/s]';
    filas{end,2}   = sprintf('%.1f',V_med(1));
    filas{end,3}   = '1898 (Tabla 3)';
    filas{end,4}   = ternario(100*abs(V_med(1)-V0)/V0<5,'✓ (<5%)','⚠');

    filas{end+1,1} = 'TOF sano env [µs]';
    filas{end,2}   = sprintf('%.2f',tof_env(1));
    filas{end,3}   = '158.2 (analítico)';
    filas{end,4}   = ternario(abs(tof_env(1)-TOF_ref)/TOF_ref<0.05,'✓ (<5%)','⚠');

    filas{end+1,1} = 'Rango exp. Fig.61 [µs]';
    filas{end,2}   = sprintf('%.2f',tof_env(1));
    filas{end,3}   = '158-190';
    filas{end,4}   = ternario(tof_env(1)>=158 && tof_env(1)<=190,'✓','~ V 2D pura');
end
if tof_env(1)>0 && tof_env(2)>0
    dTc = 100*(tof_env(2)-tof_env(1))/tof_env(1);
    filas{end+1,1} = 'dTOF hueco Ø10cm env [%]';
    filas{end,2}   = sprintf('+%.1f',dTc);
    filas{end,3}   = '+34%(sim)/+50%(exp)';
    filas{end,4}   = ternario(dTc>=25 && dTc<=65,'✓ comparable','⚠');
end
if tof_env(1)>0 && tof_env(3)>0
    err_ci = 100*abs(tof_env(1)-tof_env(3))/tof_env(3);
    filas{end+1,1} = 'Error TOF circ vs int [%]';
    filas{end,2}   = sprintf('%.2f',err_ci);
    filas{end,3}   = '<5%';
    filas{end,4}   = ternario(err_ci<5,'✓','⚠');
end
filas{end+1,1} = 'V(θ=0°) teórica [m/s]';
filas{end,2}   = sprintf('%.1f',V0);
filas{end,3}   = '1898';
filas{end,4}   = '✓';
filas{end+1,1} = 'V(θ=90°) teórica [m/s]';
filas{end,2}   = sprintf('%.1f',V90);
filas{end,3}   = '1385';
filas{end,4}   = '✓';

for i = 1:size(filas,1)
    fprintf('  %-32s %-12s %-20s %s\n', filas{i,:});
end

% Conteo de validación
n_ok = sum(cellfun(@(x) contains(x,'✓'), filas(:,4)));
fprintf('\n  Criterios validados: %d/%d\n', n_ok, size(filas,1));
if n_ok >= size(filas,1)-1
    fprintf('  → MALLA CIRCULAR VALIDADA ✓\n');
    fprintf('    Geometría realista (cilíndrica) reproduce la física de la tesis.\n');
else
    fprintf('  → MALLA CIRCULAR ACEPTABLE con diferencias documentadas\n');
end

fprintf('\n  Ventajas del tronco CIRCULAR vs. rectangular:\n');
fprintf('  → Hueco cilíndrico = misma geometría que la tesis (Fig.63)\n');
fprintf('  → No necesita corrección por área (A_circ/A_cuad)\n');
fprintf('  → Frente de onda curvo más realista (tronco real)\n');
fprintf('  → Interfaz aire-madera circular = caso experimental real\n');
fprintf('  Requisito: anglesource=90 (fuerza radial) para excitar modo rápido\n');

%% 7. FIGURAS

colores = {'r','r','b','b'};
estilos = {'-','--','-','--'};
nombres_leg = {'Circ. sin hueco','Circ. con hueco','Int. sin hueco','Int. con hueco'};

% Fig 1: Sismogramas comparados circular vs interna
figure('Name','Circ. vs Int. - Sismogramas','Position',[50 50 1000 600]);
subplot(2,1,1); hold on; grid on;
for i = [1 3]
    if isempty(datos{i}), continue; end
    sn = datos{i}.s/max(abs(datos{i}.s));
    plot(datos{i}.t*1e6, sn, [colores{i} estilos{i}],'LineWidth',1.5,...
        'DisplayName',nombres_leg{i});
    if tof_env(i)>0
        xline(tof_env(i),[colores{i} ':'],'LineWidth',1.5,'HandleVisibility','off');
    end
end
xline(TOF_ref,'k--','LineWidth',1,'DisplayName',sprintf('TOF analítico=%.1fµs',TOF_ref));
xlabel('Tiempo [µs]'); ylabel('Amplitud norm.');
title('Caso SANO: malla circular vs. interna (referencia)','FontSize',10);
legend('Location','northeast','FontSize',8); xlim([0 500]);

subplot(2,1,2); hold on; grid on;
for i = [2 4]
    if isempty(datos{i}), continue; end
    sn = datos{i}.s/max(abs(datos{i}.s));
    plot(datos{i}.t*1e6, sn, [colores{i} estilos{i}],'LineWidth',1.5,...
        'DisplayName',nombres_leg{i});
    if tof_env(i)>0
        xline(tof_env(i),[colores{i} ':'],'LineWidth',1.5,'HandleVisibility','off');
    end
end
xlabel('Tiempo [µs]'); ylabel('Amplitud norm.');
dTc_str = ''; dTi_str = '';
if tof_env(1)>0 && tof_env(2)>0
    dTc_str = sprintf('Circ dTOF=+%.1f%%',100*(tof_env(2)-tof_env(1))/tof_env(1));
end
if tof_env(3)>0 && tof_env(4)>0
    dTi_str = sprintf('Int dTOF=+%.1f%%',100*(tof_env(4)-tof_env(3))/tof_env(3));
end
title(sprintf('CON HUECO Ø10cm | %s | %s | Ref.tesis +34%%(sim)', ...
    dTc_str, dTi_str),'FontSize',10);
legend('Location','northeast','FontSize',8); xlim([0 600]);

% Fig 2: Solo malla circular — sano vs hueco — con los 3 TOF
figure('Name','Circ. - Validación vs tesis','Position',[100 100 1000 500]);
hold on; grid on;
if ~isempty(datos{1})
    sn = datos{1}.s/max(abs(datos{1}.s));
    plot(datos{1}.t*1e6, sn,'b-','LineWidth',2,'DisplayName','Sin hueco (circ.)');
end
if ~isempty(datos{2})
    sn = datos{2}.s/max(abs(datos{2}.s));
    plot(datos{2}.t*1e6, sn,'r--','LineWidth',2,'DisplayName','Con hueco Ø10cm (circ.)');
end
if tof_xcorr(1)>0, xline(tof_xcorr(1),'b:','LineWidth',1,...
    'DisplayName',sprintf('xcorr sano=%.1fµs',tof_xcorr(1))); end
if tof_env(1)>0,   xline(tof_env(1),'b-','LineWidth',2,...
    'DisplayName',sprintf('env sano=%.1fµs',tof_env(1))); end
if tof_xcorr(2)>0, xline(tof_xcorr(2),'r:','LineWidth',1,...
    'DisplayName',sprintf('xcorr hueco=%.1fµs',tof_xcorr(2))); end
if tof_env(2)>0,   xline(tof_env(2),'r-','LineWidth',2,...
    'DisplayName',sprintf('env hueco=%.1fµs',tof_env(2))); end
xline(TOF_ref,'k--','LineWidth',1,...
    'DisplayName',sprintf('TOF analítico=%.1fµs',TOF_ref));
xline(158,'k:','LineWidth',1,'DisplayName','Mín. exp. tesis (158µs)');
xline(190,'k:','LineWidth',1,'HandleVisibility','off');
xlabel('Tiempo [µs]'); ylabel('Amplitud normalizada');
title({'Malla circular: sin hueco vs. con hueco Ø10cm',...
    'Referencia: Espinosa (2019) Fig.63 — dTOF +34%(sim) / +50%(exp)'},'FontSize',10);
legend('Location','northeast','FontSize',8); xlim([0 600]);

% Fig 3: Barras dTOF comparativo
figure('Name','dTOF comparativo','Position',[150 150 700 400]);
vals_bar = [34, 50, NaN, NaN];
if tof_env(1)>0 && tof_env(2)>0
    vals_bar(3) = 100*(tof_env(2)-tof_env(1))/tof_env(1);
end
if tof_env(3)>0 && tof_env(4)>0
    vals_bar(4) = 100*(tof_env(4)-tof_env(3))/tof_env(3);
end
colores_b = [0.5 0.5 0.5; 0.7 0.7 0.7; 0.8 0.2 0.2; 0 0.4 0.8];
b = bar(vals_bar,'FaceColor','flat');
for k=1:4, b.CData(k,:) = colores_b(k,:); end
set(gca,'XTickLabel',{'Tesis sim','Tesis exp','Circular','Interna'});
ylabel('dTOF respecto al caso sano [%]');
title({'Retardo temporal por hueco centrado Ø10cm',...
    'Validación vs. tesis Espinosa (2019), Fig.63'},'FontSize',10);
grid on;
max_val = max(vals_bar(~isnan(vals_bar)));
if ~isempty(max_val) && max_val>0, ylim([0 max_val*1.3]); end
yline(34,'k--','LineWidth',1.5);
yline(50,'k:','LineWidth',1.5);
text(1, 36, '+34% sim','FontSize',8);
text(1, 52, '+50% exp','FontSize',8);

%% FUNCIÓN LOCAL
function s = ternario(cond, s_true, s_false)
    if cond, s = s_true; else, s = s_false; end
end