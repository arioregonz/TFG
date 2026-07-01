%
% MALLA CIRCULAR: SANO vs HUECO(AIRE) vs HUECO(AGUA)
%  Validación contra tesis Espinosa (2019) y malla interna
%
%  Referencia: Espinosa et al. (2019), Ultrasonics 91:242-251
%    Fig.63 (pág.79): dTOF hueco aire +34%(sim)/+50%(exp)
%    El hueco con AGUA simula pudrición húmeda — Vp_agua=1480 m/s
%    cercana a V_madera → menor contraste → dTOF intermedio
%
%  Geometría: tronco Ø30cm, hueco Ø10cm centrado
%  Fuente (0.40,0.55) anglesource=90 | Receptor (0.70,0.55) | dist=0.30m

clear; clc; close all;

fprintf('\n  MALLA CIRCULAR: SANO vs HUECO(AIRE) vs HUECO(AGUA)\n');
fprintf('  Validación vs. tesis Espinosa (2019) e interna\n\n');

%% 1. PARÁMETROS

rho = 661;        % kg/m³
c11 = 1.537e9;    % Pa  (ER pino, Tabla 2)
c33 = 8.29e8;     % Pa  (ET pino, Tabla 2)
c55 = 1.81e8;     % Pa  (GRT pino, Tabla 2)
c13 = 4.289e8;    % Pa  (νRT·√(ER·ET))

V0  = sqrt(c11/rho);    % radial 1525 m/s
V90 = sqrt(c33/rho);    % tangencial 1120 m/s
dist = 0.30;
TOF_ref = dist/V0*1e6;
DT = 2.0e-7;

% Propiedades de los rellenos del defecto
Vp_aire = 343;    % m/s
Vp_agua = 1480;   % m/s
fprintf('--- Propiedades de los medios ---\n');
fprintf('  Madera pino (radial): V=%.0f m/s\n', V0);
fprintf('  Aire:  Vp=%.0f m/s (contraste alto → fuerte reflexión)\n', Vp_aire);
fprintf('  Agua:  Vp=%.0f m/s (contraste bajo → reflexión parcial)\n', Vp_agua);
fprintf('  Z_madera/Z_aire = %.0f | Z_madera/Z_agua = %.1f\n', ...
    (rho*V0)/(1.2*Vp_aire), (rho*V0)/(1000*Vp_agua));

%% 2. ARCHIVOS

archivos = {
    'Mesh_externa/Pino_C_sin_hueco_N_2/OUTPUT_FILES/AA.S0001.BXX.semd', 'Sano (circ.)';
    'Mesh_externa/Pino_C_con_hueco_N_3/OUTPUT_FILES/AA.S0001.BXX.semd', 'Hueco aire (circ.)';
    'Mesh_externa/Pino_C_con_hueco_agua/OUTPUT_FILES/AA.S0001.BXX.semd',    'Hueco agua (circ.)';
    'Mesh_interna/Pino_sin_hueco_Norm/OUTPUT_FILES/AA.S0001.BXX.semd',  'Sano (interna ref.)';
    'Mesh_interna/Pino_con_hueco_N_C/OUTPUT_FILES/AA.S0001.BXX.semd',   'Hueco aire (interna)';
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

%% 4. DETECCIÓN TOF

vw=max(round(1/(Fc*DT)),10);
tof_xcorr=zeros(nA,1); tof_th=zeros(nA,1); tof_env=zeros(nA,1);
V_med=zeros(nA,1); energia=zeros(nA,1);

fprintf('\n--- Detección de TOF (tres métodos, Cap.3 tesis) ---\n');
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
    i0=round(0.130e-3/DT); i1=min(round(0.500e-3/DT),N);
    ie=find(env(i0:i1)>0.02*max(env),1,'first');
    if ~isempty(ie), tof_env(i)=(ie+i0-1)*DT*1e6; V_med(i)=dist/(tof_env(i)*1e-6); end
    % energía de la señal recibida
    esn_e = s/max(abs(datos{1}.s));  % normalizar al máximo del SANO
    energia(i)=trapz(t, esn_e.^2);
    fprintf('  %-22s  %8.2f  %8.2f  %8.2f  %7.1f\n', ...
        datos{i}.nombre, tof_xcorr(i), tof_th(i), tof_env(i), V_med(i));
end

%% 5. ANÁLISIS dTOF Y ENERGÍA

fprintf('\n--- dTOF y atenuación por tipo de defecto ---\n');
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
end

%% 6. VALIDACIÓN vs TESIS

fprintf('\n--- Validación vs. tesis (Espinosa 2019) ---\n');

% Sano
if tof_env(1)>0
    err_V=100*abs(V_med(1)-V0)/V0;
    fprintf('  Sano circular:\n');
    fprintf('    TOF_env=%.2fµs | V_med=%.1f m/s | error V=%.1f%% %s\n', ...
        tof_env(1), V_med(1), err_V, ternario(err_V<5,'✓','⚠'));
end

% Hueco aire (referencia directa Fig.63)
if tof_env(1)>0 && tof_env(2)>0
    dTOF_aire=100*(tof_env(2)-tof_env(1))/tof_env(1);
    fprintf('  Hueco AIRE Ø10cm:\n');
    fprintf('    dTOF=+%.1f%% | Ref.tesis Fig.63: +34%%(sim)/+50%%(exp)\n', dTOF_aire);
    fprintf('    Estado: %s\n', ternario(dTOF_aire>=25 && dTOF_aire<=65, ...
        'COMPARABLE ✓','DOCUMENTAR ⚠'));
end

% Hueco agua (interpretación física)
if tof_env(1)>0 && tof_env(3)>0
    dTOF_agua=100*(tof_env(3)-tof_env(1))/tof_env(1);
    fprintf('  Hueco AGUA Ø10cm (pudrición húmeda):\n');
    fprintf('    dTOF=+%.1f%%\n', dTOF_agua);
    fprintf('    El agua (Vp=1480 m/s) tiene menor contraste con la madera\n');
    fprintf('    → onda parcialmente transmitida → dTOF MENOR que con aire\n');
    if tof_env(2)>0
        fprintf('    Comparación: dTOF_agua=+%.1f%% < dTOF_aire=+%.1f%% %s\n', ...
            dTOF_agua, 100*(tof_env(2)-tof_env(1))/tof_env(1), ...
            ternario(dTOF_agua<100*(tof_env(2)-tof_env(1))/tof_env(1),'✓ coherente','⚠'));
    end
end

% Comparación con interna
if tof_env(1)>0 && tof_env(4)>0
    err_ci=100*abs(tof_env(1)-tof_env(4))/tof_env(4);
    fprintf('  Circular vs. interna (sano):\n');
    fprintf('    TOF_circ=%.2fµs | TOF_int=%.2fµs | error=%.1f%% %s\n', ...
        tof_env(1), tof_env(4), err_ci, ternario(err_ci<5,'✓','⚠'));
end

%% 7. TABLA RESUMEN

fprintf('\n  TABLA RESUMEN\n');
fprintf('  %-24s %-10s %-12s %-12s\n','Caso','TOF[µs]','dTOF[%]','E/E_sano[%]');
fprintf('  %s\n',repmat('-',1,60));
for i=1:3
    if tof_env(i)>0
        dt=0; if i>1 && tof_env(1)>0, dt=100*(tof_env(i)-tof_env(1))/tof_env(1); end
        er=0; if tof_env(1)>0, er=100*energia(i)/energia(1); end
        fprintf('  %-24s %-10.2f %+-12.1f %-12.1f\n', datos{i}.nombre, tof_env(i), dt, er);
    end
end

fprintf('\n  CONCLUSIÓN FÍSICA:\n');
fprintf('  → Hueco con aire: máximo contraste acústico → mayor dTOF y atenuación\n');
fprintf('  → Hueco con agua: contraste moderado → defecto detectable pero atenuado\n');
fprintf('  → Relevante para distinguir pudrición seca (aire) de húmeda (agua)\n');
fprintf('  → La tesis usa hueco con aire (Fig.63); el agua es extensión del TFG\n');

%% 8. FIGURAS

% Fig 1: Sismogramas circular - tres casos
figure('Name','Circ - Sano/Aire/Agua','Position',[50 50 1000 500]);
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
title({'Malla circular: Sano vs. Hueco(aire) vs. Hueco(agua) Ø10cm', ...
    'Ref. tesis Fig.63 (hueco aire): +34%(sim)/+50%(exp)'},'FontSize',10);
legend('Location','northeast','FontSize',8); xlim([0 600]); ylim([-1.1 1.1]);

% Fig 2: Barras dTOF y energía
figure('Name','Circ - dTOF y energía','Position',[100 100 900 400]);
subplot(1,2,1);
dtof_vals=[0 0 0];
for i=1:3, if tof_env(i)>0 && tof_env(1)>0, dtof_vals(i)=100*(tof_env(i)-tof_env(1))/tof_env(1); end; end
b1=bar(dtof_vals,'FaceColor','flat');
b1.CData=[0 0.4 0.8; 0.8 0.2 0.2; 0.2 0.7 0.3];
set(gca,'XTickLabel',{'Sano','Hueco aire','Hueco agua'});
ylabel('dTOF [%]'); title('Retardo temporal por defecto','FontSize',10);
yline(34,'k--','LineWidth',1.5); yline(50,'k:','LineWidth',1.5);
text(0.6,36,'+34% sim','FontSize',8); text(0.6,52,'+50% exp','FontSize',8);
grid on;

subplot(1,2,2);
e_vals=[100 0 0];
for i=2:3, if tof_env(1)>0, e_vals(i)=100*energia(i)/energia(1); end; end
b2=bar(e_vals,'FaceColor','flat');
b2.CData=[0 0.4 0.8; 0.8 0.2 0.2; 0.2 0.7 0.3];
set(gca,'XTickLabel',{'Sano','Hueco aire','Hueco agua'});
ylabel('Energía relativa [%]'); title('Atenuación por defecto','FontSize',10);
grid on;

% Fig 3: Comparación circular vs interna (caso aire)
figure('Name','Circ vs Int - hueco aire','Position',[150 150 900 400]);
hold on; grid on;
if ~isempty(datos{2})
    sn=datos{2}.s/max(abs(datos{2}.s));
    plot(datos{2}.t*1e6, sn,'r-','LineWidth',1.5,'DisplayName','Circular (hueco aire)');
end
if ~isempty(datos{5})
    sn=datos{5}.s/max(abs(datos{5}.s));
    plot(datos{5}.t*1e6, sn,'b--','LineWidth',1.5,'DisplayName','Interna (hueco aire)');
end
xlabel('Tiempo [µs]'); ylabel('Amplitud norm.');
title('Validación cruzada: hueco aire circular vs. interna','FontSize',10);
legend('Location','northeast','FontSize',8); xlim([0 600]);

%% FUNCIÓN LOCAL
function s = ternario(cond, s_true, s_false)
    if cond, s = s_true; else, s = s_false; end
end