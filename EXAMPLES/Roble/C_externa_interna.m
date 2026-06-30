%
%  SCRIPT 1 (v2) — MALLA CIRCULAR: SANO vs HUECO(AIRE) vs HUECO(AGUA)
%  CORRECCIÓN: cálculo de energía robusto (normalizado)
%
%  Referencia: Espinosa et al. (2019), Ultrasonics 91:242-251
%    Fig.63: dTOF hueco aire +34%(sim)/+50%(exp)
%    El agua simula pudrición húmeda — Vp_agua=1480 m/s
%
%  Geometría: tronco Ø30cm, hueco Ø10cm centrado
%  Fuente (0.40,0.55) anglesource=90 | Receptor (0.70,0.55) | dist=0.30m
%
%  NOTA SOBRE LA ENERGÍA:
%    Las amplitudes absolutas de cada simulación difieren por el factor
%    de fuente y la geometría, por lo que NO se pueden comparar energías
%    brutas. Se mide la energía de cada señal NORMALIZADA A SU PROPIO
%    MÁXIMO, en una ventana fija tras la primera llegada. Así la energía
%    refleja la DURACIÓN/DISPERSIÓN del paquete recibido (cuánto se
%    "ensancha" la señal por reverberación y difracción), no su amplitud.

clear; clc; close all;

fprintf('\n  MALLA CIRCULAR: SANO vs HUECO(AIRE) vs HUECO(AGUA)\n');
fprintf('  Validación vs. tesis Espinosa (2019) e interna\n\n');

%% 1. PARÁMETROS (ROBLE — cambiar para pino)

rho = 706; c11 = 2.54e9; c33 = 1.35e9; c55 = 3.19e8; c13 = 7.5832e8;
V0  = sqrt(c11/rho);    % radial
V90 = sqrt(c33/rho);    % tangencial
dist = 0.30;
TOF_ref = dist/V0*1e6;
DT = 2.0e-7;

Vp_aire = 343; Vp_agua = 1480;
fprintf('--- Propiedades de los medios ---\n');
fprintf('  Madera (radial): V=%.0f m/s\n', V0);
fprintf('  Aire:  Vp=%.0f m/s (contraste alto → fuerte reflexión)\n', Vp_aire);
fprintf('  Agua:  Vp=%.0f m/s (contraste bajo → reflexión parcial)\n', Vp_agua);
fprintf('  Z_madera/Z_aire = %.0f | Z_madera/Z_agua = %.1f\n', ...
    (rho*V0)/(1.2*Vp_aire), (rho*V0)/(1000*Vp_agua));

%% 2. ARCHIVOS

archivos = {
    'Mesh_externa/Roble_C_sin_hueco_N_2/OUTPUT_FILES/AA.S0001.BXX.semd', 'Sano (circ.)';
    'Mesh_externa/Roble_C_con_hueco_N_3/OUTPUT_FILES/AA.S0001.BXX.semd', 'Hueco aire (circ.)';
    'Mesh_externa/Roble_C_con_hueco_agua/OUTPUT_FILES/AA.S0001.BXX.semd',    'Hueco agua (circ.)';
    'Mesh_interna/Roble_sin_hueco_Norm/OUTPUT_FILES/AA.S0001.BXX.semd',  'Sano (interna ref.)';
    'Mesh_interna/Roble_con_hueco_N_C/OUTPUT_FILES/AA.S0001.BXX.semd',   'Hueco aire (interna)';
};

fprintf('\n--- Cargando sismogramas ---\n');
nA = size(archivos,1);
datos = cell(nA,1);
for i=1:nA
    if exist(archivos{i,1},'file')
        A=load(archivos{i,1});
        datos{i}.t=A(:,1)-A(1,1); datos{i}.s=A(:,2); datos{i}.nombre=archivos{i,2};
        fprintf('  ✓ %s (%d pts)\n', archivos{i,2}, length(A));
    else
        fprintf('  ✗ NO EXISTE: %s\n', archivos{i,1});
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

%% 4. DETECCIÓN TOF Y ENERGÍA (CORREGIDA)

vw=max(round(1/(Fc*DT)),10);
tof_xcorr=zeros(nA,1); tof_th=zeros(nA,1); tof_env=zeros(nA,1);
V_med=zeros(nA,1); energia=zeros(nA,1);

% Ventana fija para medir energía (tras la primera llegada posible)
i0_e = round(0.130e-3/DT);
i1_e = min(round(0.600e-3/DT), N);

fprintf('\n--- Detección de TOF (tres métodos) y energía ---\n');
fprintf('  %-22s  xcorr[µs] thresh[µs]  env[µs]  V[m/s]\n','Caso');
fprintf('  %s\n',repmat('-',1,66));

for i=1:nA
    if isempty(datos{i}), continue; end
    t=datos{i}.t; s=datos{i}.s; sn=s/max(abs(s));
    % xcorr
    [c,lags]=xcorr(sn,chirp_full,'coeff'); [~,ic]=max(c);
    tof_xcorr(i)=lags(ic)*DT*1e6;
    % threshold
    sig=std(sn(1:round(0.100e-3/DT))); umb=max(8*sig,0.10);
    ith=find(abs(sn)>umb,1,'first'); if ~isempty(ith), tof_th(i)=t(ith)*1e6; end
    % envolvente
    env=sqrt(movmean(sn.^2,vw));
    ie=find(env(i0_e:i1_e)>0.02*max(env),1,'first');
    if ~isempty(ie), tof_env(i)=(ie+i0_e-1)*DT*1e6; V_med(i)=dist/(tof_env(i)*1e-6); end
    % ENERGÍA CORREGIDA: sobre señal normalizada a SU PROPIO máximo,
    % en ventana fija → mide la dispersión/duración del paquete
    energia(i)=trapz(t(i0_e:i1_e), sn(i0_e:i1_e).^2);
    fprintf('  %-22s  %8.2f  %8.2f  %8.2f  %7.1f\n', ...
        datos{i}.nombre, tof_xcorr(i), tof_th(i), tof_env(i), V_med(i));
end

%% 5. ANÁLISIS dTOF Y ENERGÍA

fprintf('\n--- dTOF y energía relativa por tipo de defecto ---\n');
fprintf('  %-22s  dTOF[%%]   E/E_sano[%%]  Interpretación\n','Caso');
fprintf('  %s\n',repmat('-',1,70));

if tof_env(1)>0
    E_sano = energia(1);
    casos = {2,'Hueco aire','reflexión total (Z alto)'; ...
             3,'Hueco agua','reflexión parcial (Z medio)'};
    for k=1:size(casos,1)
        idx=casos{k,1};
        if idx<=nA && tof_env(idx)>0
            dtof=100*(tof_env(idx)-tof_env(1))/tof_env(1);
            erel=100*energia(idx)/E_sano;
            fprintf('  %-22s  %+7.1f   %8.1f    %s\n', ...
                casos{k,2}, dtof, erel, casos{k,3});
        end
    end
    fprintf('  (energía sobre señal normalizada → mide dispersión del paquete)\n');
end

%% 6. VALIDACIÓN vs TESIS

fprintf('\n--- Validación vs. tesis (Espinosa 2019) ---\n');
if tof_env(1)>0
    err_V=100*abs(V_med(1)-V0)/V0;
    fprintf('  Sano circular: TOF=%.2fµs | V=%.1f m/s | error=%.1f%% %s\n', ...
        tof_env(1), V_med(1), err_V, ternario(err_V<5,'✓','⚠'));
end
if tof_env(1)>0 && tof_env(2)>0
    dTOF_aire=100*(tof_env(2)-tof_env(1))/tof_env(1);
    fprintf('  Hueco AIRE: dTOF=+%.1f%% | tesis +34%%(sim)/+50%%(exp) %s\n', ...
        dTOF_aire, ternario(dTOF_aire>=25 && dTOF_aire<=65,'✓','⚠'));
end
if tof_env(1)>0 && tof_env(3)>0
    dTOF_agua=100*(tof_env(3)-tof_env(1))/tof_env(1);
    dTOF_aire=100*(tof_env(2)-tof_env(1))/tof_env(1);
    fprintf('  Hueco AGUA: dTOF=+%.1f%% (< aire +%.1f%%) %s\n', ...
        dTOF_agua, dTOF_aire, ternario(dTOF_agua<dTOF_aire,'✓ coherente','⚠'));
    fprintf('    El agua (Z≈Z_madera) transmite parcialmente → menor retardo\n');
end

%% 7. TABLA RESUMEN

fprintf('\n  TABLA RESUMEN\n');
fprintf('  %-22s %-10s %-10s %-12s\n','Caso','TOF[µs]','dTOF[%]','E/E_sano[%]');
fprintf('  %s\n',repmat('-',1,56));
for i=1:3
    if tof_env(i)>0
        dt=0; if i>1, dt=100*(tof_env(i)-tof_env(1))/tof_env(1); end
        er=100*energia(i)/energia(1);
        fprintf('  %-22s %-10.2f %+-10.1f %-12.1f\n', datos{i}.nombre, tof_env(i), dt, er);
    end
end

fprintf('\n  CONCLUSIÓN FÍSICA:\n');
fprintf('  → Aire: máximo contraste → mayor dTOF\n');
fprintf('  → Agua: contraste moderado → defecto detectable pero atenuado\n');
fprintf('  → Distingue pudrición seca (aire) de húmeda (agua)\n');

%% 8. FIGURAS

% Fig 1: Sismogramas
figure('Name','Circ - Sano/Aire/Agua','Position',[50 50 1000 500],'Color','w');
hold on; grid on;
cc = {'b','r','g'};
for i=1:3
    if isempty(datos{i}), continue; end
    sn=datos{i}.s/max(abs(datos{i}.s));
    plot(datos{i}.t*1e6, sn, cc{i}, 'LineWidth',1.5, 'DisplayName',datos{i}.nombre);
    if tof_env(i)>0, xline(tof_env(i),[cc{i} ':'],'LineWidth',1.5,'HandleVisibility','off'); end
end
xline(TOF_ref,'k--','LineWidth',1,'DisplayName',sprintf('TOF analítico=%.1fµs',TOF_ref));
xlabel('Tiempo [µs]'); ylabel('Amplitud norm.');
title('Malla circular: Sano vs. Hueco(aire) vs. Hueco(agua) Ø10cm','FontSize',10);
legend('Location','northeast','FontSize',8); xlim([0 600]); ylim([-1.1 1.1]);

% Fig 2: Barras dTOF y energía (CORREGIDA)
figure('Name','Circ - dTOF y energía','Position',[100 100 900 400],'Color','w');

subplot(1,2,1);
dtof_vals=[0 0 0];
for i=1:3, if tof_env(i)>0 && tof_env(1)>0, dtof_vals(i)=100*(tof_env(i)-tof_env(1))/tof_env(1); end; end
b1=bar(dtof_vals,'FaceColor','flat');
b1.CData=[0 0.4 0.8; 0.8 0.2 0.2; 0.2 0.7 0.3];
set(gca,'XTickLabel',{'Sano','Hueco aire','Hueco agua'});
ylabel('dTOF [%]'); title('Retardo temporal por defecto','FontSize',10);
yline(34,'k--','LineWidth',1.5);
text(0.55,36,'+34% sim','FontSize',8);
grid on; ylim([0 max(dtof_vals)*1.25]);

subplot(1,2,2);
e_vals=[100 0 0];
for i=2:3, if energia(1)>0, e_vals(i)=100*energia(i)/energia(1); end; end
b2=bar(e_vals,'FaceColor','flat');
b2.CData=[0 0.4 0.8; 0.8 0.2 0.2; 0.2 0.7 0.3];
set(gca,'XTickLabel',{'Sano','Hueco aire','Hueco agua'});
ylabel('Energía relativa [%]');
title('Energía del paquete recibido','FontSize',10);
grid on; ylim([0 max(e_vals)*1.25]);

%% FUNCIÓN LOCAL
function s = ternario(cond, s_true, s_false)
    if cond, s = s_true; else, s = s_false; end
end