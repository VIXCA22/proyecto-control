clear; clc; close all;

%% === 1. Cargar tablas conservando nombres originales ===
fertig   = readtable('3-4 fertig.csv',  'PreserveVariableNames', true);
chien    = readtable('3-4 chen.csv',    'PreserveVariableNames', true);
sintesis = readtable('3-4 sintesis analitica.csv', 'PreserveVariableNames', true);

%% === 2. Renombrar columnas al estilo "Pandas" para las 3 tablas ===
fertig = renamevars(fertig, ...
    {'Input (rpm)', 'Setpoint (rpm)', 'Output (PWM)'}, ...
    {'Referencia (rpm)', 'Señal de Control (PWM)', 'Salida (rpm)'});

chien = renamevars(chien, ...
    {'Input (rpm)', 'Setpoint (rpm)', 'Output (PWM)'}, ...
    {'Referencia (rpm)', 'Señal de Control (PWM)', 'Salida (rpm)'});

sintesis = renamevars(sintesis, ...
    {'Input (rpm)', 'Setpoint (rpm)', 'Output (PWM)'}, ...
    {'Referencia (rpm)', 'Señal de Control (PWM)', 'Salida (rpm)'});


%% =========== FERTIG: NORMALIZACIÓN ESPECIAL (190/160/220 -> 120/90/120) ===========

[t_f, ref_f, sc_f, y_f] = normalizar_fertig(fertig);

figure; hold on; grid on; box on;
plot(t_f, ref_f, 'k--', 'LineWidth', 1.6);                          % referencia
plot(t_f, sc_f,  '--',  'Color',[0.3 0.3 0.3], 'LineWidth',1.4);    % señal de control
plot(t_f, y_f,   'b',   'LineWidth', 1.6);                          % salida
xlabel('Tiempo (s)');
ylabel('Señal (rpm)');
title('Fertig — Señales normalizadas (120 \rightarrow 90 \rightarrow 120 rpm)');
legend('Referencia', 'Señal de Control', 'Salida');


%% =========== CHIEN: reescalado a referencia 90–120 (LINEAL) ===========

[t_c, ref_c, sc_c, y_c] = reescalar_a_90_120(chien);

figure; hold on; grid on; box on;
plot(t_c, ref_c, 'k--', 'LineWidth', 1.6);
plot(t_c, sc_c, '--', 'Color',[0.3 0.3 0.3], 'LineWidth',1.4);
plot(t_c, y_c, 'r', 'LineWidth', 1.6);
xlabel('Tiempo (s)');
ylabel('Señal (rpm)');
title('Chien — Señales reescaladas a 90–120 rpm');
legend('Referencia', 'Señal de Control', 'Salida');


%% ======= SÍNTESIS ANALÍTICA: recentrado a referencia 90–120 (OFFSET) =======
[t_s, ref_s, sc_s, y_s] = recentrar_a_90_120(sintesis);

figure; hold on; grid on; box on;
plot(t_s, ref_s, 'k--', 'LineWidth', 1.6);
plot(t_s, sc_s, '--', 'Color',[0.3 0.3 0.3], 'LineWidth',1.4);
plot(t_s, y_s, 'm', 'LineWidth', 1.6);
xlabel('Tiempo (s)');
ylabel('Señal (rpm)');
title('Síntesis Analítica — Señales re-centradas a 90–120 rpm');
legend('Referencia', 'Señal de Control', 'Salida');


%% ================== FUNCIÓN: NORMALIZAR FERTIG (MAPPING POR NIVELES) ==================
% Fertig tiene niveles 190, 160, 220 en la referencia original.
% Queremos:
%   190 -> 120
%   160 -> 90
%   220 -> 120
% y que la salida se mueva aproximadamente igual.

function [t_norm, ref_new, sc_new, y_new] = normalizar_fertig(T)

    t   = T.("Tiempo (s)");
    ref = T.("Referencia (rpm)");
    sc  = T.("Señal de Control (PWM)");
    y   = T.("Salida (rpm)");

    % Copia de referencia
    ref_new = ref;

    % Mapping por niveles (ajusta aquí si ves otros valores)
    mask190 = (ref == 190);
    mask160 = (ref == 160);
    mask220 = (ref == 220);

    % 190 -> 120
    ref_new(mask190) = 120;
    % 160 -> 90
    ref_new(mask160) = 90;
    % 220 -> 120
    ref_new(mask220) = 120;

    % Desplazar SC y Y con el mismo delta que su referencia
    delta = ref_new - ref;      % cuánto se movió cada punto de referencia

    sc_new = sc + delta;
    y_new  = y  + delta;

    % Normalizar tiempo al primer cambio en la referencia nueva
    idx_step = find(ref_new ~= ref_new(1), 1, 'first');
    if isempty(idx_step)
        t0 = t(1);
    else
        t0 = t(idx_step);
    end

    t_norm = t - t0;
end


%% ================== FUNCIÓN: recentrar (OFFSET) ==================
function [t_norm, ref_new, sc_new, y_new] = recentrar_a_90_120(T)
    t   = T.("Tiempo (s)");
    ref = T.("Referencia (rpm)");
    sc  = T.("Señal de Control (PWM)");
    y   = T.("Salida (rpm)");

    idx_step = find(abs(ref - ref(1)) > 0.5, 1);
    if isempty(idx_step)
        t0 = t(1); 
        ref_base = ref(1);
    else
        t0 = t(idx_step);
        ref_base = mean(ref(1:idx_step-1));
    end

    t_norm = t - t0;
    offset = 90 - ref_base;

    ref_new = ref + offset;
    sc_new  = sc  + offset;
    y_new   = y   + offset;
end


%% ================== FUNCIÓN: reescalar (LINEAL) ==================
function [t_norm, ref_new, sc_new, y_new] = reescalar_a_90_120(T)
    t   = T.("Tiempo (s)");
    ref = T.("Referencia (rpm)");
    sc  = T.("Señal de Control (PWM)");
    y   = T.("Salida (rpm)");

    r_min = min(ref);
    r_max = max(ref);

    a = 30 / (r_max - r_min);
    b = 90 - a * r_min;

    ref_new = a * ref + b;
    sc_new  = a * sc  + b;
    y_new   = a * y   + b;

    r_mid = (min(ref_new) + max(ref_new)) / 2;
    idx_step = find(ref_new > r_mid, 1);

    if isempty(idx_step)
        t0 = t(1);
    else
        t0 = t(idx_step);
    end

    t_norm = t - t0;
end

%% === DEFINIR SEÑALES DE CONTROL ===
u_fertig   = sc_f;    
u_chien    = sc_c;    
u_sintesis = sc_s;    


%% === CALCULAR TVu Y Umax ===

% FERTIG
TVud_fertig   = sum(abs(diff(u_fertig)));
Umaxd_fertig  = max(abs(u_fertig));

% CHIEN
TVud_chien    = sum(abs(diff(u_chien)));
Umaxd_chien   = max(abs(u_chien));

% SÍNTESIS
TVud_sintesis = sum(abs(diff(u_sintesis)));
Umaxd_sintesis = max(abs(u_sintesis));


%% === ÍNDICES REGULATORIOS USANDO VALOR EN RÉGIMEN PREVIO ===
fprintf('\n================ ÍNDICES REGULATORIOS (definición estándar) ================\n');
fprintf('%-18s | %10s | %10s | %10s | %10s\n', ...
        'Controlador','IAEd','ISEd','Tuvd','Emax');
fprintf('-----------------------------------------------------------------------------\n');

%% ---------- FERTIK: ventana regulatoria 25–30 s ----------
% Ventana de análisis (regulatorio)
idxF = (t_f >= 25) & (t_f <= 30);
tF   = t_f(idxF);
yF   = y_f(idxF);
uF   = sc_f(idxF);

% Estimar y_ss justo antes de la perturbación, por ejemplo 24–25 s
idxF_ss = (t_f >= 24) & (t_f < 25);
yF_ss   = mean(y_f(idxF_ss));   % valor deseado en regulatorio

eF = yF - yF_ss;                % error regulatorio

IAEd_fertig = trapz(tF, abs(eF));
ISEd_fertig = trapz(tF, eF.^2);
Tuvd_fertig = sum(abs(diff(uF)));
Emax_fertig = max(abs(eF));

fprintf('%-18s | %10.4f | %10.4f | %10.4f | %10.4f\n', ...
        'Fertik y Sharpe', IAEd_fertig, ISEd_fertig, Tuvd_fertig, Emax_fertig);

%% ---------- CHIEN: ventana regulatoria 58–62 s ----------
idxC = (t_c >= 58) & (t_c <= 62);
tC   = t_c(idxC);
yC   = y_c(idxC);
uC   = sc_c(idxC);

% Régimen antes de la perturbación, por ejemplo 57–58 s
idxC_ss = (t_c >= 57) & (t_c < 58);
yC_ss   = mean(y_c(idxC_ss));

eC = yC - yC_ss;

IAEd_chien = trapz(tC, abs(eC));
ISEd_chien = trapz(tC, eC.^2);
Tuvd_chien = sum(abs(diff(uC)));
Emax_chien = max(abs(eC));

fprintf('%-18s | %10.4f | %10.4f | %10.4f | %10.4f\n', ...
        'Chien (1988)', IAEd_chien, ISEd_chien, Tuvd_chien, Emax_chien);

%% ---------- SÍNTESIS: ventana regulatoria 28–32 s ----------
idxS = (t_s >= 28) & (t_s <= 32);
tS   = t_s(idxS);
yS   = y_s(idxS);
uS   = sc_s(idxS);

% Régimen antes de la perturbación, por ejemplo 27–28 s
idxS_ss = (t_s >= 27) & (t_s < 28);
yS_ss   = mean(y_s(idxS_ss));

eS = yS - yS_ss;

IAEd_sintesis  = trapz(tS, abs(eS));
ISEd_sintesis  = trapz(tS, eS.^2);
Tuvd_sintesis  = sum(abs(diff(uS)));
Emax_sintesis  = max(abs(eS));

fprintf('%-18s | %10.4f | %10.4f | %10.4f | %10.4f\n', ...
        'Síntesis Analítica', IAEd_sintesis, ISEd_sintesis, Tuvd_sintesis, Emax_sintesis);

fprintf('-----------------------------------------------------------------------------\n');
