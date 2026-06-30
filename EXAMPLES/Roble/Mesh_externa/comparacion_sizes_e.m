%
%  MALLA CIRCULAR: BARRIDO DE TAMAÑOS DE HUECO
%  Hueco de aire centrado, diámetros Ø10/Ø5/Ø2 cm
%  Validación vs. tesis Espinosa (2019) y malla interna
%
%  Referencia: Espinosa et al. (2019)
%    Fig.42 (pág.62): dTOF aumenta con el tamaño del defecto
%    Fig.63 (pág.79): hueco Ø10cm → +34%(sim)/+50%(exp)
%
%  Geometría: tronco Ø30cm, hueco centrado de diámetro variable
%  Fuente (0.40,0.55) anglesource=90 | Receptor (0.70,0.55) | dist=0.30m

clear; clc; close all;

fprintf('\n  MALLA CIRCULAR: BARRIDO DE TAMAÑOS DE HUECO\n');
fprintf('  Hueco aire centrado Ø10/Ø5/Ø2 cm | vs. tesis e interna\n\n');

%% 1. PARÁMETROS

rho=706; c11=2.54e9; c33=1.35e9;
V0=sqrt(c11/rho); dist=0.30; TOF_ref=dist/V0*1e6; DT=2.0e-7;

%% 2. ARCHIVOS — orden: sano, Ø10, Ø5, Ø2 (circular) + interna referencia

archivos = {
    'Mesh_externa/Roble_C_sin_hueco_N_2/OUTPUT_FILES/AA.S0001.BXX.semd', 'Sano',        0.0;
    'Mesh_externa/Roble_C_con_hueco_N_3/OUTPUT_FILES/AA.S0001.BXX.semd', 'Hueco Ø10cm', 10.0;
    'Mesh_externa/Roble_C_hueco_5/OUTPUT_FILES/AA.S0001.BXX.semd',     'Hueco Ø5cm',  5.0;
    'Mesh_externa/Roble_C_hueco_2/OUTPUT_FILES/AA.S0001.BXX.semd',     'Hueco Ø2cm',  2.0;
};

% Referencia interna (hueco cuadrado, para comparar tendencia)
archivos_int = {
    'Mesh_interna/Roble_con_hueco_N_C/OUTPUT_FILES/AA.S0001.BXX.semd',          '10x10cm', 10.0;
    'Mesh_interna/Roble_con_hueco_N_C_Mediano/OUTPUT_FILES/AA.S0001.BXX.semd',  '5x5cm',   5.0;
    'Mesh_interna/Roble_con_hueco_N_C_Pequeño/OUTPUT_FILES/AA.S0001.BXX.semd',  '2x2cm',   2.0;
};

fprintf('--- Cargando sismogramas (malla circular) ---\n');
nC=size(archivos,1);
datos=cell(nC,1);
for i=1:nC
    if exist(archivos{i,1},'file')
        A=load(archivos{i,1});
        datos{i}.t=A(:,1)-A(1,1); datos{i}.s=A(:,2);
        datos{i}.nombre=archivos{i,2}; datos{i}.diam=archivos{i,3};
        fprintf('  ✓ %s (%d pts)\n', archivos{i,2}, length(A));
    else
        fprintf('  ✗ NO EXISTE: %s\n', archivos{i,1});
        datos{i}=[];
    end
end

fprintf('--- Cargando referencia interna (hueco cuadrado) ---\n');
nI=size(archivos_int,1);
datos_int=cell(nI,1);
for i=1:nI
    if exist(archivos_int{i,1},'file')
        A=load(archivos_int{i,1});
        datos_int{i}.t=A(:,1)-A(1,1); datos_int{i}.s=A(:,2);
        datos_int{i}.nombre=archivos_int{i,2}; datos_int{i}.diam=archivos_int{i,3};
        fprintf('  ✓ %s (%d pts)\n', archivos_int{i,2}, length(A));
    else
        fprintf('  ✗ NO EXISTE: %s\n', archivos_int{i,1});
        datos_int{i}=[];
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

%% 4. DETECCIÓN TOF (circular)

vw=max(round(1/(Fc*DT)),10);
tof_env=zeros(nC,1); tof_xc=zeros(nC,1); tof_th=zeros(nC,1);
diam=zeros(nC,1);

fprintf('\n--- TOF malla circular (tres métodos) ---\n');
fprintf('  %-14s  Ø[cm]  xcorr[µs] thresh[µs] env[µs]\n','Caso');
fprintf('  %s\n',repmat('-',1,58));

for i=1:nC
    if isempty(datos{i}), continue; end
    t=datos{i}.t; sn=datos{i}.s/max(abs(datos{i}.s)); diam(i)=datos{i}.diam;
    [c,lags]=xcorr(sn,chirp_full,'coeff'); [~,ic]=max(c); tof_xc(i)=lags(ic)*DT*1e6;
    sig=std(sn(1:round(0.100e-3/DT))); umb=max(8*sig,0.10);
    ith=find(abs(sn)>umb,1,'first'); if ~isempty(ith), tof_th(i)=t(ith)*1e6; end
    env=sqrt(movmean(sn.^2,vw));
    i0=round(0.130e-3/DT); i1=min(round(0.500e-3/DT),N);
    ie=find(env(i0:i1)>0.02*max(env),1,'first');
    if ~isempty(ie), tof_env(i)=(ie+i0-1)*DT*1e6; end
    fprintf('  %-14s  %5.1f  %8.2f  %8.2f  %7.2f\n', ...
        datos{i}.nombre, diam(i), tof_xc(i), tof_th(i), tof_env(i));
end

%% 5. TOF INTERNA (referencia tendencia)

tof_env_int=zeros(nI,1); lado_int=zeros(nI,1);
% Necesita el TOF sano interno
A_si=load('Mesh_interna/Roble_sin_hueco_Norm/OUTPUT_FILES/AA.S0001.BXX.semd');
s_si=A_si(:,2)/max(abs(A_si(:,2)));
env_si=sqrt(movmean(s_si.^2,vw));
i0=round(0.130e-3/DT); i1=min(round(0.500e-3/DT),length(s_si));
ie=find(env_si(i0:i1)>0.02*max(env_si),1,'first');
tof_sano_int=(ie+i0-1)*DT*1e6;

for i=1:nI
    if isempty(datos_int{i}), continue; end
    sn=datos_int{i}.s/max(abs(datos_int{i}.s)); lado_int(i)=datos_int{i}.diam;
    env=sqrt(movmean(sn.^2,vw));
    ie=find(env(i0:min(i1,length(sn)))>0.02*max(env),1,'first');
    if ~isempty(ie), tof_env_int(i)=(ie+i0-1)*DT*1e6; end
end

%% 6. dTOF vs TAMAÑO

fprintf('\n--- dTOF vs. tamaño del defecto ---\n');
fprintf('  CIRCULAR (hueco cilíndrico, vs. tesis Fig.63):\n');
fprintf('  %-14s  Ø[cm]  TOF[µs]  dTOF[%%]  Ref.tesis\n','Caso');
fprintf('  %s\n',repmat('-',1,56));
ref_t={'base','+34%(sim)','~+22%(sim)','~+11%(sim)'};
for i=1:nC
    if tof_env(i)>0
        dt=0; if i>1 && tof_env(1)>0, dt=100*(tof_env(i)-tof_env(1))/tof_env(1); end
        fprintf('  %-14s  %5.1f  %7.2f  %+6.1f   %s\n', ...
            datos{i}.nombre, diam(i), tof_env(i), dt, ref_t{i});
    end
end

fprintf('\n  INTERNA (hueco cuadrado, referencia tendencia):\n');
fprintf('  %-14s lado[cm] TOF[µs]  dTOF[%%]\n','Caso');
fprintf('  %s\n',repmat('-',1,48));
fprintf('  %-14s  %5.1f  %7.2f  %+6.1f\n','Sano',0,tof_sano_int,0);
for i=1:nI
    if tof_env_int(i)>0
        dt=100*(tof_env_int(i)-tof_sano_int)/tof_sano_int;
        fprintf('  %-14s  %5.1f  %7.2f  %+6.1f\n', ...
            datos_int{i}.nombre, lado_int(i), tof_env_int(i), dt);
    end
end

%% 7. VERIFICACIÓN TENDENCIA

fprintf('\n--- Verificación de tendencia (Fig.42 tesis) ---\n');
dtof_circ=zeros(nC,1);
for i=1:nC, if tof_env(i)>0 && tof_env(1)>0, dtof_circ(i)=100*(tof_env(i)-tof_env(1))/tof_env(1); end; end
monotona = all(diff(dtof_circ(diam>0 & tof_env>0)) ...
    .*sign(diff(diam(diam>0 & tof_env>0))) >= -5);
fprintf('  Tendencia esperada: dTOF aumenta con Ø del defecto\n');
fprintf('  dTOF(Ø10)=%.1f%% > dTOF(Ø5)=%.1f%% > dTOF(Ø2)=%.1f%%\n', ...
    dtof_circ(2), dtof_circ(3), dtof_circ(4));
if dtof_circ(2)>=dtof_circ(3) && dtof_circ(3)>=dtof_circ(4)
    fprintf('  → TENDENCIA CORRECTA ✓ (monótona creciente con tamaño)\n');
else
    fprintf('  → Revisar: tendencia no perfectamente monótona ⚠\n');
end

%% 8. FIGURAS

% Fig 1: Sismogramas circular por tamaño
figure('Name','Circ - Barrido tamaños','Position',[50 50 1000 500]);
hold on; grid on;
cc={'k','r','b','g'};
for i=1:nC
    if isempty(datos{i}), continue; end
    sn=datos{i}.s/max(abs(datos{i}.s));
    plot(datos{i}.t*1e6, sn, cc{i},'LineWidth',1.5,'DisplayName',datos{i}.nombre);
    if tof_env(i)>0, xline(tof_env(i),[cc{i} ':'],'LineWidth',1,'HandleVisibility','off'); end
end
xline(TOF_ref,'k--','LineWidth',1,'DisplayName',sprintf('TOF analítico=%.1fµs',TOF_ref));
xlabel('Tiempo [µs]'); ylabel('Amplitud norm.');
title('Malla circular: barrido de tamaños de hueco (aire, centrado)','FontSize',10);
legend('Location','northeast','FontSize',8); xlim([0 600]);

% Fig 2: dTOF vs diámetro — circular vs interna
figure('Name','dTOF vs tamaño','Position',[100 100 800 500]);
hold on; grid on;
d_circ=diam(tof_env>0 & diam>0); dt_circ=dtof_circ(tof_env>0 & diam>0);
[d_circ,is]=sort(d_circ); dt_circ=dt_circ(is);
plot(d_circ, dt_circ,'ro-','LineWidth',2,'MarkerSize',8,...
    'MarkerFaceColor','r','DisplayName','Circular (hueco cilíndrico)');

d_int=lado_int(tof_env_int>0); dt_int=zeros(sum(tof_env_int>0),1);
cnt=0;
for i=1:nI
    if tof_env_int(i)>0
        cnt=cnt+1;
        dt_int(cnt)=100*(tof_env_int(i)-tof_sano_int)/tof_sano_int;
    end
end
[d_int,is]=sort(d_int); dt_int=dt_int(is);
plot(d_int, dt_int,'bs--','LineWidth',2,'MarkerSize',8,...
    'MarkerFaceColor','b','DisplayName','Interna (hueco cuadrado)');

% Referencias tesis
plot(10, 34,'k^','MarkerSize',12,'MarkerFaceColor','k',...
    'DisplayName','Tesis Ø10cm +34% (sim)');
plot(10, 50,'kv','MarkerSize',12,'MarkerFaceColor','y',...
    'DisplayName','Tesis Ø10cm +50% (exp)');
xlabel('Diámetro / lado del defecto [cm]'); ylabel('dTOF [%]');
title({'dTOF vs. tamaño del defecto','Validación tendencia Fig.42/63 tesis'},'FontSize',10);
legend('Location','northwest','FontSize',8);
xlim([0 12]);

%% FUNCIÓN LOCAL
function s = ternario(cond, s_true, s_false)
    if cond, s = s_true; else, s = s_false; end
end