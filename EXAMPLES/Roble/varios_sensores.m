%
%  ANÁLISIS DE LOCALIZACIÓN CON 16 SENSORES — Roble circular
%  Tomografía simplificada: patrón angular de TOF
%
%  Referencia: Espinosa et al. (2019), Cap.4
%    16 sensores alrededor del tronco (pág.56)
%    El patrón angular de dTOF localiza el defecto
%
%  Geometría: tronco Ø30cm, centro (0.55,0.55)
%    16 sensores en circunferencia (radio 0.145m, cada 22.5°)
%    Fuente en ángulo 180° (x=0.40, borde izquierdo), anglesource=90
%    S01=0° (enfrentado a fuente) ... S09=180° (junto a fuente)
%
%  Casos:
%    1. Sano (referencia)
%    2. Hueco centrado Ø10cm → retardo simétrico en sensores opuestos
%    3. Hueco excéntrico → retardo localizado en sensores tras el hueco

clear; clc; close all;

fprintf('\n  ANÁLISIS DE LOCALIZACIÓN CON 16 SENSORES — Roble\n');
fprintf('  Patrón angular de TOF para localizar el defecto\n\n');

%% 1. PARÁMETROS Y GEOMETRÍA

rho=706; c11=2.54e9; c33=1.35e9; c55=3.19e8; c13=7.5832e8;
V0=sqrt(c11/rho); V90=sqrt(c33/rho);
DT=2.0e-7;
cx=0.55; cz=0.55; R=0.15; R_sensor=0.145;
n_sens=16;

% Posición de la fuente (ángulo 180°, borde izquierdo)
xs_fuente = 0.40; zs_fuente = 0.55;

% Posiciones de los 16 sensores
ang_sens = (0:n_sens-1) * 360/n_sens;  % grados
xs_sens = cx + R_sensor*cosd(ang_sens);
zs_sens = cz + R_sensor*sind(ang_sens);

% Posición de los huecos
hueco_centro = [0.55, 0.55];   % centrado
hueco_exc    = [0.50, 0.60];   % excéntrico (UL, arriba-izquierda)

fprintf(' Configuración \n');
fprintf('  Tronco Ø30cm, centro (%.2f,%.2f)\n', cx, cz);
fprintf('  Fuente en (%.2f,%.2f) ángulo 180°\n', xs_fuente, zs_fuente);
fprintf('  16 sensores cada 22.5° en radio %.3fm\n', R_sensor);
fprintf('  Hueco centrado: (%.2f,%.2f)\n', hueco_centro(1), hueco_centro(2));
fprintf('  Hueco excéntrico: (%.2f,%.2f)\n', hueco_exc(1), hueco_exc(2));

%% 2. CARGA DE LOS 3 CASOS

casos = {
    'Mesh_externa/Roble_C_16sensores_sano/OUTPUT_FILES',   'Sano';
    'Mesh_externa/Roble_C_16sensores_hueco/OUTPUT_FILES',  'Hueco centrado';
    'Mesh_externa/Roble_C_16sensores_exc/OUTPUT_FILES',    'Hueco excéntrico';
};
n_casos = size(casos,1);

% Matriz de TOF: filas=sensores, columnas=casos
TOF = nan(n_sens, n_casos);

% Chirp para xcorr
Ts=45e-6; f0c=22e3; f1c=50e3; Fc=36e3;
Nc=round(Ts/DT); t_c=(0:Nc-1)*DT;
f_t=(f1c-f0c)*t_c/Ts+f0c; phase=2*pi*f_t.*t_c;
win=exp(-((t_c-Ts/2).^2)/(2*(Ts/6)^2));
chirp_s=cos(phase).*win;
vw=max(round(1/(Fc*DT)),10);

fprintf('\nCargando sismogramas de cada caso \n');
for c = 1:n_casos
    carpeta = casos{c,1};
    n_ok = 0;
    for s = 1:n_sens
        archivo = sprintf('%s/AA.S%02d.BXX.semd', carpeta, s);
        if exist(archivo,'file')
            A = load(archivo);
            t = A(:,1)-A(1,1);
            sig = A(:,2);
            N = length(t);
            sn = sig/max(abs(sig));
            % TOF por envolvente
            env = sqrt(movmean(sn.^2, vw));
            i0 = round(0.130e-3/DT); i1 = min(round(0.500e-3/DT),N);
            ie = find(env(i0:i1)>0.02*max(env),1,'first');
            if ~isempty(ie)
                TOF(s,c) = (ie+i0-1)*DT*1e6;
            end
            n_ok = n_ok + 1;
        end
    end
    fprintf('  %s: %d/%d sensores cargados\n', casos{c,2}, n_ok, n_sens);
end

%% 3. dTOF POR SENSOR (respecto al sano)

fprintf('\n dTOF por sensor (respecto al caso sano) \n');
fprintf('  Sensor  Ángulo   TOF_sano  TOF_centro  dTOF_c   TOF_exc  dTOF_e\n');
fprintf('  %s\n', repmat('-',1,68));

dTOF_centro = nan(n_sens,1);
dTOF_exc    = nan(n_sens,1);

for s = 1:n_sens
    if ~isnan(TOF(s,1))
        if ~isnan(TOF(s,2))
            dTOF_centro(s) = 100*(TOF(s,2)-TOF(s,1))/TOF(s,1);
        end
        if ~isnan(TOF(s,3))
            dTOF_exc(s) = 100*(TOF(s,3)-TOF(s,1))/TOF(s,1);
        end
        fprintf('  S%02d     %5.1f°   %7.2f   %7.2f    %+5.1f%%   %7.2f  %+5.1f%%\n', ...
            s, ang_sens(s), TOF(s,1), TOF(s,2), dTOF_centro(s), ...
            TOF(s,3), dTOF_exc(s));
    end
end

%% 4. LOCALIZACIÓN DEL DEFECTO

fprintf('\nLocalización del defecto\n');

% Hueco centrado: máximo dTOF debe estar en el sensor enfrentado a la fuente
[~, s_max_c] = max(dTOF_centro);
fprintf('  Hueco CENTRADO:\n');
fprintf('    Máximo dTOF en sensor S%02d (ángulo %.1f°) = %+.1f%%\n', ...
    s_max_c, ang_sens(s_max_c), dTOF_centro(s_max_c));
fprintf('    Esperado: S01 (0°, enfrentado a fuente) — rayo cruza el centro\n');

% Hueco excéntrico: máximo dTOF localiza la posición angular del defecto
[~, s_max_e] = max(dTOF_exc);
fprintf('  Hueco EXCÉNTRICO (0.50,0.60):\n');
fprintf('    Máximo dTOF en sensor S%02d (ángulo %.1f°) = %+.1f%%\n', ...
    s_max_e, ang_sens(s_max_e), dTOF_exc(s_max_e));
% Calcular el ángulo esperado: línea fuente→hueco prolongada
ang_hueco = atan2d(hueco_exc(2)-zs_fuente, hueco_exc(1)-xs_fuente);
if ang_hueco<0, ang_hueco=ang_hueco+360; end
fprintf('    Dirección fuente→hueco: %.1f° (el defecto desvía estos rayos)\n', ang_hueco);

%% 5. FIGURA 1: Patrón polar de dTOF

figure('Name','Patrón angular dTOF','Position',[50 50 1100 550]);

subplot(1,2,1);
ang_rad = deg2rad([ang_sens, ang_sens(1)]);  % cerrar el círculo
dt_c = [dTOF_centro; dTOF_centro(1)]';
dt_c(isnan(dt_c)) = 0;
polarplot(ang_rad, dt_c, 'r-o', 'LineWidth',2, 'MarkerFaceColor','r');
title({'Hueco CENTRADO','dTOF simétrico (máx. en sensor enfrentado)'},'FontSize',10);
ax=gca; ax.ThetaZeroLocation='right'; ax.ThetaDir='counterclockwise';

subplot(1,2,2);
dt_e = [dTOF_exc; dTOF_exc(1)]';
dt_e(isnan(dt_e)) = 0;
polarplot(ang_rad, dt_e, 'b-o', 'LineWidth',2, 'MarkerFaceColor','b');
title({'Hueco EXCÉNTRICO','dTOF localizado (máx. en dirección del defecto)'},'FontSize',10);
ax=gca; ax.ThetaZeroLocation='right'; ax.ThetaDir='counterclockwise';

%% 6. FIGURA 2: Mapa espacial del tronco con sensores

figure('Name','Mapa tronco con 16 sensores','Position',[100 100 700 700]);
hold on; axis equal; grid on;

% Tronco
th=linspace(0,2*pi,200);
plot(cx+R*cos(th), cz+R*sin(th),'k-','LineWidth',2,'HandleVisibility','off');

% Sensores coloreados por dTOF del caso excéntrico
for s=1:n_sens
    dt = dTOF_exc(s);
    if isnan(dt), dt=0; end
    col = [min(dt/max(dTOF_exc(~isnan(dTOF_exc))),1) 0 ...
           max(1-dt/max(dTOF_exc(~isnan(dTOF_exc))),0)];
    plot(xs_sens(s), zs_sens(s),'o','MarkerSize',12,...
        'MarkerFaceColor',col,'MarkerEdgeColor','k','HandleVisibility','off');
    text(xs_sens(s)+0.01*cosd(ang_sens(s)), zs_sens(s)+0.01*sind(ang_sens(s)), ...
        sprintf('S%02d',s),'FontSize',7);
end

% Fuente
plot(xs_fuente, zs_fuente,'g^','MarkerSize',16,'MarkerFaceColor','g',...
    'DisplayName','Fuente');

% Hueco excéntrico
plot(hueco_exc(1)+0.05*cos(th), hueco_exc(2)+0.05*sin(th),'b--','LineWidth',2,...
    'DisplayName','Hueco excéntrico');

% Líneas fuente→sensores con mayor dTOF
for s=1:n_sens
    if ~isnan(dTOF_exc(s)) && dTOF_exc(s) > 0.5*max(dTOF_exc(~isnan(dTOF_exc)))
        plot([xs_fuente xs_sens(s)],[zs_fuente zs_sens(s)],'r:','LineWidth',1,...
            'HandleVisibility','off');
    end
end

xlabel('x [m]'); ylabel('z [m]');
title({'Localización del defecto excéntrico','Rojo=mayor retardo (rayo cruza el hueco)'},...
    'FontSize',10);
legend('Location','northeastoutside','FontSize',9);
xlim([0.35 0.75]); ylim([0.35 0.75]);

%% 7. FIGURA 3: dTOF vs ángulo (cartesiano)

figure('Name','dTOF vs ángulo','Position',[150 150 900 450]);
hold on; grid on;
plot(ang_sens, dTOF_centro,'r-o','LineWidth',2,'MarkerFaceColor','r',...
    'DisplayName','Hueco centrado');
plot(ang_sens, dTOF_exc,'b-s','LineWidth',2,'MarkerFaceColor','b',...
    'DisplayName','Hueco excéntrico');
xline(0,'k:','DisplayName','Sensor enfrentado (0°)');
xline(180,'g:','DisplayName','Junto a fuente (180°)');
xlabel('Ángulo del sensor [°]'); ylabel('dTOF [%]');
title('Patrón angular de dTOF: centrado vs excéntrico','FontSize',10);
legend('Location','northeast','FontSize',9);
xlim([0 360]); xticks(0:45:360);

%% FUNCIÓN LOCAL
function s = ternario(cond, s_true, s_false)
    if cond, s = s_true; else, s = s_false; end
end