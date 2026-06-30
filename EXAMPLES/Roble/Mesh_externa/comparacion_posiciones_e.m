%
%   MALLA CIRCULAR: BARRIDO DE POSICIONES DE HUECO
%  Hueco de aire Ø10cm en 5 posiciones: Centro, UR, UL, DL, DR
%  Validación vs. tesis Espinosa (2019) y malla interna
%
%  Referencia: Espinosa et al. (2019)
%    Cap.4, pág.63-84: dependencia del dTOF con la posición del defecto
%    pág.84: dTOF_centro > dTOF_excéntrico (defecto en la ruta del rayo)
%    Fig.44-51: mapas de dTOF según posición
%
%  Geometría: tronco Ø30cm (centro 0.55,0.55), hueco Ø10cm
%    Centro: (0.55,0.55) — intercepta el rayo directo
%    UR: (0.60,0.60) | UL: (0.50,0.60) | DL: (0.50,0.50) | DR: (0.60,0.50)
%    (offset 0.05m del centro → solo difracción, no interceptan rayo)
%  Fuente (0.40,0.55) anglesource=90 | Receptor (0.70,0.55) | dist=0.30m
%  Rayo directo: línea horizontal z=0.55

clear; clc; close all;

fprintf('\n  MALLA CIRCULAR: BARRIDO DE POSICIONES DE HUECO (Ø10cm)\n');
fprintf('  vs. tesis Espinosa (2019) Cap.4 e interna\n\n');

%% 1. PARÁMETROS

rho=706; c11=2.54e9; c33=1.35e9;
V0=sqrt(c11/rho); dist=0.30; TOF_ref=dist/V0*1e6; DT=2.0e-7;
cx=0.55; cz=0.55; R_h=0.05; z_rayo=0.55;

%% 2. ARCHIVOS — circular

% Posiciones y sus coordenadas de centro del hueco
pos_info = {
    'Sano',   NaN,  NaN;
    'Centro', 0.55, 0.55;
    'UR',     0.60, 0.60;
    'UL',     0.50, 0.60;
    'DL',     0.50, 0.50;
    'DR',     0.60, 0.50;
};

archivos = {
    'Roble_C_sin_hueco_N_2/OUTPUT_FILES/AA.S0001.BXX.semd';
    'Roble_C_con_hueco_N_3/OUTPUT_FILES/AA.S0001.BXX.semd';
    'Roble_C_hueco_UR/OUTPUT_FILES/AA.S0001.BXX.semd';
    'Roble_C_hueco_UL/OUTPUT_FILES/AA.S0001.BXX.semd';
    'Roble_C_hueco_DL/OUTPUT_FILES/AA.S0001.BXX.semd';
    'Roble_C_hueco_DR/OUTPUT_FILES/AA.S0001.BXX.semd';
};

fprintf('--- Cargando sismogramas (malla circular) ---\n');
nP=size(archivos,1);
datos=cell(nP,1);
for i=1:nP
    if exist(archivos{i},'file')
        A=load(archivos{i});
        datos{i}.t=A(:,1)-A(1,1); datos{i}.s=A(:,2);
        datos{i}.nombre=pos_info{i,1};
        datos{i}.hx=pos_info{i,2}; datos{i}.hz=pos_info{i,3};
        fprintf('  ✓ %s (%d pts)\n', pos_info{i,1}, length(A));
    else
        fprintf('  ✗ NO EXISTE: %s\n', archivos{i});
        datos{i}=[];
    end
end

%% 3. CHIRP

Ts=45e-6; f0c=22e3; f1c=50e3; Fc=36e3;
Nc=round(Ts/DT); t_c=(0:Nc-1)*DT;
f_t=(f1c-f0c)*t_c/Ts+f0c; phase=2*pi*f_t.*t_c;
win=exp(-((t_c-Ts/2).^2)/(2*(Ts/6)^2));
chirp_s=cos(phase).*win;
N=3000; if ~isempty(datos{1}), N=length(datos{1}.t); end
chirp_full=zeros(N,1); chirp_full(1:Nc)=chirp_s;

%% 4. DETECCIÓN TOF

vw=max(round(1/(Fc*DT)),10);
tof_env=zeros(nP,1); tof_xc=zeros(nP,1); tof_th=zeros(nP,1);
intercepta=false(nP,1);

fprintf('\n--- Detección de TOF (tres métodos) ---\n');
fprintf('  %-8s  xcorr[µs] thresh[µs]  env[µs]  Rayo?\n','Pos');
fprintf('  %s\n',repmat('-',1,52));

for i=1:nP
    if isempty(datos{i}), continue; end
    t=datos{i}.t; sn=datos{i}.s/max(abs(datos{i}.s));
    [c,lags]=xcorr(sn,chirp_full,'coeff'); [~,ic]=max(c); tof_xc(i)=lags(ic)*DT*1e6;
    sig=std(sn(1:round(0.100e-3/DT))); umb=max(8*sig,0.10);
    ith=find(abs(sn)>umb,1,'first'); if ~isempty(ith), tof_th(i)=t(ith)*1e6; end
    env=sqrt(movmean(sn.^2,vw));
    i0=round(0.130e-3/DT); i1=min(round(0.500e-3/DT),N);
    ie=find(env(i0:i1)>0.02*max(env),1,'first');
    if ~isempty(ie), tof_env(i)=(ie+i0-1)*DT*1e6; end
    % ¿intercepta el rayo?
    if ~isnan(datos{i}.hz)
        intercepta(i) = abs(datos{i}.hz - z_rayo) < R_h*0.90;
    end
    rayo_str = '—';
    if i>1, rayo_str = ternario(intercepta(i),'SÍ','no'); end
    fprintf('  %-8s  %8.2f  %8.2f  %7.2f   %s\n', ...
        datos{i}.nombre, tof_xc(i), tof_th(i), tof_env(i), rayo_str);
end

%% 5. dTOF POR POSICIÓN

fprintf('\n--- dTOF por posición (envolvente) ---\n');
fprintf('  %-8s  Centro hueco    TOF[µs]  dTOF[%%]  Rayo?  Ref.tesis\n','Pos');
fprintf('  %s\n',repmat('-',1,70));

dtof_pos=zeros(nP,1);
for i=1:nP
    if tof_env(i)==0, continue; end
    if i>1 && tof_env(1)>0
        dtof_pos(i)=100*(tof_env(i)-tof_env(1))/tof_env(1);
    end
    if i==1
        fprintf('  %-8s  base            %7.2f  %+6.1f   —      base\n', ...
            datos{i}.nombre, tof_env(i), 0);
    else
        rayo_str=ternario(intercepta(i),'SÍ ','no ');
        ref_str=ternario(intercepta(i),'+34%(sim,en ruta)','<5%(difracción)');
        fprintf('  %-8s  (%.2f,%.2f)     %7.2f  %+6.1f   %s   %s\n', ...
            datos{i}.nombre, datos{i}.hx, datos{i}.hz, ...
            tof_env(i), dtof_pos(i), rayo_str, ref_str);
    end
end

%% 6. INTERPRETACIÓN (tendencias tesis pág.84)

fprintf('\n--- Interpretación (Cap.4, pág.63-84 tesis) ---\n');

% Tendencia 1: centro > excéntrico
dtof_centro = dtof_pos(2);
dtof_exc = dtof_pos(3:6);
dtof_exc_max = max(dtof_exc(dtof_exc~=0));
if isempty(dtof_exc_max), dtof_exc_max=0; end
fprintf('  Tendencia 1 (pág.84): dTOF_centro > dTOF_excéntrico\n');
fprintf('    dTOF_centro=%.1f%%  dTOF_excéntrico_max=%.1f%%  → %s\n', ...
    dtof_centro, dtof_exc_max, ...
    ternario(dtof_centro>dtof_exc_max,'✓ CUMPLE','✗ revisar'));

fprintf('  Tendencia 2: solo el hueco centrado intercepta el rayo directo\n');
fprintf('    → mayor dTOF en centro (efecto directo)\n');
fprintf('    → posiciones UR/UL/DL/DR solo difractan (offset 0.05m del rayo)\n');

% Tendencia 3: asimetría por anisotropía
fprintf('  Tendencia 3: asimetría entre posiciones (anisotropía roble)\n');
if dtof_pos(3)>0 && dtof_pos(5)>0
    fprintf('    dTOF_UR=%.1f%% vs dTOF_DL=%.1f%%\n', dtof_pos(3), dtof_pos(5));
end
if dtof_pos(4)>0 && dtof_pos(6)>0
    fprintf('    dTOF_UL=%.1f%% vs dTOF_DR=%.1f%%\n', dtof_pos(4), dtof_pos(6));
end
fprintf('    La asimetría refleja V_radial>V_tangencial del roble\n');

%% 7. TABLA RESUMEN

fprintf('\n  TABLA RESUMEN — POSICIONES\n');
fprintf('  %-8s %-14s %-10s %-10s %s\n','Pos','Centro','TOF[µs]','dTOF[%]','Estado');
fprintf('  %s\n',repmat('-',1,58));
for i=1:nP
    if tof_env(i)==0, continue; end
    if i==1
        est='base';
        cstr='—';
    else
        cstr=sprintf('(%.2f,%.2f)',datos{i}.hx,datos{i}.hz);
        if intercepta(i)
            est=ternario(dtof_pos(i)>20,'✓ en ruta','⚠');
        else
            est=ternario(dtof_pos(i)<10,'✓ difracción','⚠');
        end
    end
    fprintf('  %-8s %-14s %-10.2f %+-10.1f %s\n', ...
        datos{i}.nombre, cstr, tof_env(i), dtof_pos(i), est);
end

%% 8. FIGURAS

% Fig 1: Sismogramas por posición
figure('Name','Circ - Posiciones','Position',[50 50 1000 550]);
hold on; grid on;
cc={'k','r','b','g','m','c'};
for i=1:nP
    if isempty(datos{i}), continue; end
    sn=datos{i}.s/max(abs(datos{i}.s));
    plot(datos{i}.t*1e6, sn, cc{i},'LineWidth',1.3,'DisplayName',datos{i}.nombre);
    if tof_env(i)>0, xline(tof_env(i),[cc{i} ':'],'LineWidth',1,'HandleVisibility','off'); end
end
xline(TOF_ref,'k--','LineWidth',1,'DisplayName',sprintf('TOF analítico=%.1fµs',TOF_ref));
xlabel('Tiempo [µs]'); ylabel('Amplitud norm.');
title('Malla circular: barrido de posiciones del hueco Ø10cm','FontSize',10);
legend('Location','northeast','FontSize',8); xlim([0 600]);

% Fig 2: Mapa espacial de dTOF (vista del tronco)
figure('Name','Mapa dTOF posiciones','Position',[100 100 650 600]);
hold on; axis equal; grid on;
% Tronco
th=linspace(0,2*pi,100);
plot(cx+0.15*cos(th), cz+0.15*sin(th),'k-','LineWidth',2,'HandleVisibility','off');
% Rayo directo
plot([0.40 0.70],[0.55 0.55],'b--','LineWidth',1.5,'DisplayName','Rayo directo');
plot(0.40,0.55,'g^','MarkerSize',12,'MarkerFaceColor','g','DisplayName','Fuente');
plot(0.70,0.55,'rv','MarkerSize',12,'MarkerFaceColor','r','DisplayName','Receptor');
% Huecos con color por dTOF
for i=2:nP
    if isempty(datos{i}) || isnan(datos{i}.hx), continue; end
    % Color según dTOF (rojo=alto, azul=bajo)
    dt=dtof_pos(i);
    plot(datos{i}.hx+R_h*cos(th), datos{i}.hz+R_h*sin(th),'-',...
        'Color',[min(dt/40,1) 0 max(1-dt/40,0)],'LineWidth',2,...
        'HandleVisibility','off');
    text(datos{i}.hx, datos{i}.hz, sprintf('%s\n%.1f%%',datos{i}.nombre,dt),...
        'HorizontalAlignment','center','FontSize',8,'FontWeight','bold');
end
xlabel('x [m]'); ylabel('z [m]');
title({'Mapa de dTOF por posición del hueco','Rojo=mayor retardo, Azul=menor'},'FontSize',10);
legend('Location','northeastoutside','FontSize',8);
xlim([0.35 0.75]); ylim([0.35 0.75]);

% Fig 3: Barras dTOF por posición
figure('Name','dTOF por posición','Position',[150 150 700 400]);
nombres_b={}; vals_b=[];
for i=2:nP
    if tof_env(i)>0
        nombres_b{end+1}=datos{i}.nombre;
        vals_b(end+1)=dtof_pos(i);
    end
end
b=bar(vals_b,'FaceColor','flat');
for k=1:length(vals_b)
    if vals_b(k)>20, b.CData(k,:)=[0.8 0.2 0.2];
    else, b.CData(k,:)=[0 0.4 0.8]; end
end
set(gca,'XTickLabel',nombres_b);
ylabel('dTOF [%]');
title({'dTOF por posición del defecto','Centro (en ruta) >> excéntrico (difracción)'},'FontSize',10);
yline(34,'k--','LineWidth',1.5,'DisplayName','Ref. tesis +34% (centro)');
grid on; legend('FontSize',8,'Location','northeast');

%% FUNCIÓN LOCAL
function s = ternario(cond, s_true, s_false)
    if cond, s = s_true; else, s = s_false; end
end