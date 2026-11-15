% 1. Cargar datos
datos = readmatrix('DatosMotor.csv');
tiempo = datos(:, 1);
senal_control = datos(:, 2);
salida_velocidad = datos(:, 3);

% 2. Encontrar puntos de normalización
idx_escalon = find(senal_control ~= 90, 1, 'first');
yinicial = mean(salida_velocidad(1 : idx_escalon - 1));
tiempo_escalon = tiempo(idx_escalon);

% 3. Normalizar TODAS las señales 
tiempo_norm_total = tiempo - tiempo_escalon; % t=0 será el escalón
salida_norm_total = salida_velocidad - yinicial; % Señal centrada en 0
entrada_norm_total = senal_control - 90; % Señal de 0 a 30

% 4. Aplicar Filtrado
windowSize = 6; 
salida_filt_total = movmean(salida_norm_total, windowSize);

% 5. Definir la ventana de visualización
lookback_points = 30; 
idx_start = idx_escalon - lookback_points;
if idx_start < 1
    idx_start = 1; 
end
idx_end = length(tiempo);

% 6. Crear los vectores finales para el PLOT
tiempo_plot = tiempo_norm_total(idx_start : idx_end);
salida_orig_plot = salida_norm_total(idx_start : idx_end);
salida_filt_plot = salida_filt_total(idx_start : idx_end);
entrada_plot = entrada_norm_total(idx_start : idx_end);

% 7. Visualizar


% Subplot 1: Curva de Reacción (Salida)
 % (2 filas, 1 columna, gráfico 1)
plot(tiempo_plot, salida_orig_plot, 'b-', 'DisplayName', 'Señal Normalizada');
hold on;
plot(tiempo_plot, salida_filt_plot, 'r-', 'LineWidth', 2, 'DisplayName', 'Señal Filtrada');
xline(0, 'k--', 'DisplayName', 't = 0 (Escalón)');
yline(0, 'k-');
hold off;
grid on;
ylabel('Incremento de Velocidad');
title('Curva de Reacción (Salida Normalizada y Filtrada)');
legend('Location', 'southeast');



%systemIdentification
load('P1D.mat','-mat');


%% === Bloque adicional: Identificación 123c (Alfaro) y sobreposición en el 1er subplot ===
% Usa variables ya existentes: tiempo_plot, salida_filt_plot, entrada_plot

% 1) Baselines y amplitudes (antes y después del escalón)
mask_pre  = tiempo_plot < 0;
mask_tail = true(size(tiempo_plot));
tailN = max(25, round(0.05*numel(tiempo_plot)));
if tailN < numel(tiempo_plot)
    mask_tail(1:end-tailN) = false;  % toma la "cola" final
end

y0  = mean(salida_filt_plot(mask_pre),  'omitnan');   if isnan(y0),  y0  = salida_filt_plot(1); end
u0  = mean(entrada_plot(mask_pre),      'omitnan');   if isnan(u0),  u0  = 0; end
yss = mean(salida_filt_plot(mask_tail), 'omitnan');   if isnan(yss), yss = salida_filt_plot(end); end
uss = mean(entrada_plot(mask_tail),     'omitnan');   if isnan(uss), uss = entrada_plot(end); end

dy = yss - y0;          % cambio total en salida
du = uss - u0;          % paso real en la entrada
if du == 0
    warning('Δu = 0: no se puede estimar la ganancia. Verifica los datos de entrada.');
    du = eps;
end

% 2) Niveles 25% y 75% y tiempos de cruce (interpolación lineal)
y25 = y0 + 0.25*dy;
y75 = y0 + 0.75*dy;

% --- t25 ---
idx = find(salida_filt_plot >= y25, 1, 'first');
if isempty(idx) || idx == 1
    t25 = NaN;
else
    t1 = tiempo_plot(idx-1); t2 = tiempo_plot(idx);
    y1 = salida_filt_plot(idx-1); y2 = salida_filt_plot(idx);
    t25 = t1 + (y25 - y1) * (t2 - t1) / (y2 - y1);
end

% --- t75 ---
idx = find(salida_filt_plot >= y75, 1, 'first');
if isempty(idx) || idx == 1
    t75 = NaN;
else
    t1 = tiempo_plot(idx-1); t2 = tiempo_plot(idx);
    y1 = salida_filt_plot(idx-1); y2 = salida_filt_plot(idx);
    t75 = t1 + (y75 - y1) * (t2 - t1) / (y2 - y1);
end

if any(isnan([t25, t75]))
    warning('No se pudieron determinar t25/t75; revisa la señal o aumenta el suavizado.');
end

% 3) Parámetros 123c (según la lámina)
tau = 0.9102 * (t75 - t25);
tm  = 1.2620 * t25 - 0.2620 * t75;
Kp  = dy / du;

% Saneamiento básico
if ~isfinite(tau) || tau <= 0, tau = abs(tau); warning('tau no válido; se usó |tau|.'); end
if ~isfinite(tm)  || tm  < 0, tm  = 0;         warning('tm<0; se fija a 0 (sin retardo).'); end

% 4) Modelo FOPDT y simulación sobre malla uniforme (lsim exige tiempo limpio)
s_tf = tf('s');
G = Kp/(tau*s_tf + 1);
if tm > 0, G.InputDelay = tm; end

tt = tiempo_plot(:);
uu = entrada_plot(:);
mask = isfinite(tt) & isfinite(uu);
tt = tt(mask); uu = uu(mask);

[tt, ia] = unique(tt, 'stable');  % tiempos estrictamente crecientes
uu = uu(ia);

dt = median(diff(tt));
if ~isfinite(dt) || dt <= 0
    dt = mean(diff(tt));
end
t_u = (tt(1):dt:tt(end)).';                 % malla uniforme
eu  = interp1(tt, uu, t_u, 'linear', 'extrap');

y_model = lsim(G, eu, t_u);                 % respuesta del modelo



% === (nuevo) 4.5) Calcular FIT del modelo 123c vs salida filtrada ===
yu = interp1(tiempo_plot(:), salida_filt_plot(:), t_u, 'linear','extrap');  % salida real en malla uniforme
fit_123c = 100*(1 - norm(yu - y_model)/norm(yu - mean(yu)));

% === (reemplazo) 5) Sobreponer en el primer subplot con FIT en la leyenda ===

hold on;
plot(t_u, y_model, 'k--', 'LineWidth', 1.8, ...
     'DisplayName', sprintf('Modelo 123c — Fit %.2f%%', fit_123c));
legend('show','Location','southeast');
hold off;

% (opcional) imprimir en consola
fprintf('\n=== Modelo 123c (Alfaro) ===\n');
fprintf('Fit respecto a salida filtrada: %.2f %%\n', fit_123c);


%% === Bloque adicional: Mínimos Cuadrados (ARX 1,1,1) + fit y gráfico ===
% Usa: tiempo_plot, entrada_plot, salida_filt_plot
% Modelo: y(k) = a1*y(k-1) + b1*u(k-1)

% 1) Preparar vectores columna y limpiar NaN/Inf
t = tiempo_plot(:);
u = entrada_plot(:);
y = salida_filt_plot(:);
mask = isfinite(t) & isfinite(u) & isfinite(y);
t = t(mask); u = u(mask); y = y(mask);
N = numel(y);

% 2) Construir matrices para LS (ARX(na=1, nb=1, nk=1))
nk = 1;                 % retraso puro de 1 muestra
k0 = 1 + nk;            % primer índice válido para ecuación
X  = [ y(k0-1:N-1) , u(k0-nk:N-nk) ];   % [y(k-1)  u(k-1)]
Y  = y(k0:N);                             

% 3) Resolver por mínimos cuadrados (θ = X\Y)
theta = X\Y;            % equivalente a (X'*X)\(X'*Y) pero más estable
a1 = theta(1);
b1 = theta(2);

% 4) Predicción recursiva sin lsim (evita problemas de malla no uniforme)
y_LS = zeros(N,1);
y_LS(1) = y(1);
for k = k0:N
    y_LS(k) = a1*y_LS(k-1) + b1*u(k-nk);
end

% 5) Fit (%) respecto a la salida filtrada
fit_LS = 100*(1 - norm(y - y_LS)/norm(y - mean(y)));

% 6) Mostrar el modelo discreto equivalente y el fit
fprintf('\n=== Mínimos Cuadrados (ARX 1,1,1) ===\n');
fprintf('a1 = %.6f, b1 = %.6f, nk = %d\n', a1, b1, nk);
fprintf('G(z) = (%.6f) z^{-%d} / (1 - %.6f z^{-1})\n', b1, nk, a1);
try
    % Mostrar como tf (solo para impresión; no se usa lsim aquí)
    Ts_est = median(diff(t(isfinite(t))));
    if ~isfinite(Ts_est) || Ts_est <= 0, Ts_est = 0; end
    num = [zeros(1,nk) b1];
    den = [1 -a1];
    Gz_LS = tf(num, den, Ts_est);
    disp(Gz_LS);
catch
    % si no hay Control System Toolbox, se omite
end
fprintf('Fit (respecto a salida filtrada): %.2f %%\n', fit_LS);

% 7) Añadir al PRIMER subplot

hold on;
plot(t, y_LS, 'm-.', 'LineWidth', 1.6, ...
    'DisplayName', sprintf('Mínimos Cuadrados (ARX 1,1,1) — Fit %.2f%%', fit_LS));
legend('show','Location','southeast');
hold off;


%% === Señal y FIT del P1D importado + sobreponer en 1er subplot ===
% Requiere en el workspace: P1D, tiempo_plot, entrada_plot, salida_filt_plot

if exist('P1D','var')
    % 1) Convertir a tf continuo (si es idproc/idtf)
    if isa(P1D,'idproc')
        G_P1D_tf = tf(idtf(P1D));
    elseif isa(P1D,'idtf')
        G_P1D_tf = tf(P1D);
    elseif isa(P1D,'tf')
        G_P1D_tf = P1D;
    else
        error('El modelo P1D no es de tipo idproc, idtf o tf.');
    end

    % 2) Asegurar malla de tiempo UNIFORME para simular/compare
    tt = tiempo_plot(:);
    uu = entrada_plot(:);
    yy = salida_filt_plot(:);
    m = isfinite(tt) & isfinite(uu) & isfinite(yy);
    tt = tt(m); uu = uu(m); yy = yy(m);
    [tt, ia] = unique(tt,'stable'); uu = uu(ia); yy = yy(ia);
    dt = median(diff(tt)); if ~isfinite(dt) || dt <= 0, dt = mean(diff(tt)); end
    t_u = (tt(1):dt:tt(end)).';
    eu  = interp1(tt, uu, t_u, 'linear','extrap');
    yu  = interp1(tt, yy, t_u, 'linear','extrap');

    % 3) FIT con el método del Toolbox (preferido) y señal simulada
    fit_txt = 'N/A';
    try
        if license('test','System_Identification_Toolbox')
            data_id = iddata(yu, eu, dt);
            [yhat, fitInfo] = compare(data_id, G_P1D_tf);   % respuesta y % oficial
            % Extraer salida simulada desde iddata devuelta
            if isprop(yhat,'OutputData')
                y_p1d = yhat.OutputData; 
            else
                y_p1d = yhat.y;  % compatibilidad
            end
            y_p1d = y_p1d(:);
            fit_val = fitInfo(1).Fit.Percentage;           
            fit_txt = sprintf('%.2f', fit_val);
        else
            error('No Toolbox'); % forzar vía manual
        end
    catch
        % Ruta manual si compare no está disponible
        y_p1d = lsim(G_P1D_tf, eu, t_u);
        fit_val = 100*(1 - norm(yu - y_p1d)/norm(yu - mean(yu)));
        fit_txt = sprintf('%.2f', fit_val);
    end

    % 4) Graficar sobre el PRIMER subplot con el FIT en la leyenda

    hold on;
    plot(t_u, y_p1d, '-', 'LineWidth', 1.8, ...
        'DisplayName', sprintf('P1D (Toolbox) — Fit %s%%', fit_txt));
    legend('show','Location','southeast');
    hold off;

    % 5) Reporte en consola
    fprintf('\n=== Señal y FIT del P1D importado ===\n');
    fprintf('FIT (%%): %s\n', fit_txt);
else
    warning('No se encontró el modelo P1D en el workspace.');
end



% 6) Imprimir resultados
fprintf('\n=== Método 123c (Alfaro) ===\n');
fprintf('Kp = %.6g\n', Kp);
fprintf('tau = %.6g s\n', tau);
fprintf('tm  = %.6g s\n', tm);
fprintf('Δu  = %.6g,  Δy = %.6g\n', du, dy);
fprintf('G(s) = %.6g / (%.6g s + 1) * e^{-%g s}\n', Kp, tau, tm);
disp(G);

%% === Bloque adicional: Método de Ho et al. (35%–85%) ===
% Requiere: tiempo_plot, salida_filt_plot, entrada_plot (ya creadas)

% 1) Baseline y amplitudes (antes y después del escalón)
mask_pre  = tiempo_plot < 0;
mask_tail = true(size(tiempo_plot));
tailN = max(25, round(0.05*numel(tiempo_plot)));
if tailN < numel(tiempo_plot), mask_tail(1:end-tailN) = false; end

y0  = mean(salida_filt_plot(mask_pre),  'omitnan'); if isnan(y0),  y0  = salida_filt_plot(1); end
u0  = mean(entrada_plot(mask_pre),      'omitnan'); if isnan(u0),  u0  = 0; end
yss = mean(salida_filt_plot(mask_tail), 'omitnan'); if isnan(yss), yss = salida_filt_plot(end); end
uss = mean(entrada_plot(mask_tail),     'omitnan'); if isnan(uss), uss = entrada_plot(end); end

dy = yss - y0;                    % cambio total de salida
du = uss - u0;                    % amplitud real del paso
if du == 0, warning('Δu = 0; se fuerza eps.'); du = eps; end
Kp = dy/du;

% 2) Tiempos t35 y t85 por interpolación lineal
pct1 = 0.35; pct2 = 0.85;
y35  = y0 + pct1*dy;
y85  = y0 + pct2*dy;

% t35
idx = find(salida_filt_plot >= y35, 1, 'first');
if isempty(idx) || idx==1
    t35 = NaN;
else
    t1 = tiempo_plot(idx-1); t2 = tiempo_plot(idx);
    y1 = salida_filt_plot(idx-1); y2 = salida_filt_plot(idx);
    t35 = t1 + (y35 - y1)*(t2 - t1)/(y2 - y1);
end
% t85
idx = find(salida_filt_plot >= y85, 1, 'first');
if isempty(idx) || idx==1
    t85 = NaN;
else
    t1 = tiempo_plot(idx-1); t2 = tiempo_plot(idx);
    y1 = salida_filt_plot(idx-1); y2 = salida_filt_plot(idx);
    t85 = t1 + (y85 - y1)*(t2 - t1)/(y2 - y1);
end

if any(isnan([t35,t85]))
    warning('No se pudieron determinar t35/t85; revisa señal o filtrado.');
end

% 3) Parámetros de Ho et al. (a=0.670, b=1.290)
a = 0.670;  b = 1.290;
T  = a*(t85 - t35);
L  = b*t35 + (1-b)*t85;        % = 1.29*t35 - 0.29*t85

% Saneamiento básico
if ~isfinite(T) || T<=0, T = abs(T); warning('T no válido; se usó |T|.'); end
if ~isfinite(L) || L<0,  L = max(0,L); warning('L<0; se satura a 0.'); end

% 4) Modelo FOPDT de Ho et al. y simulación (malla uniforme)
s_tf = tf('s');
G_ho = Kp/(T*s_tf + 1);
if L > 0, G_ho.InputDelay = L; end

% malla uniforme (si no existe de bloques previos, la creamos)
if ~(exist('t_u','var') && exist('eu','var'))
    tt = tiempo_plot(:); uu = entrada_plot(:);
    m = isfinite(tt) & isfinite(uu);
    tt = tt(m); uu = uu(m);
    [tt, ia] = unique(tt,'stable'); uu = uu(ia);
    dt = median(diff(tt)); if ~isfinite(dt) || dt<=0, dt = mean(diff(tt)); end
    t_u = (tt(1):dt:tt(end)).';
    eu  = interp1(tt, uu, t_u, 'linear','extrap');
end
y_ho = lsim(G_ho, eu, t_u);

% 5) FIT respecto a la salida filtrada y leyenda
yu = interp1(tiempo_plot(:), salida_filt_plot(:), t_u, 'linear','extrap');
fit_ho = 100*(1 - norm(yu - y_ho)/norm(yu - mean(yu)));



% 7) Impresión en consola
fprintf('\n=== Método Ho et al. (35%%–85%%) ===\n');
fprintf('Kp = %.6g, T = %.6g s, L = %.6g s\n', Kp, T, L);
fprintf('G(s) = %.6g / (%.6g s + 1) * e^{-%g s}\n', Kp, T, L);
disp(G_ho);

%% === Bloque final: IAE, ISE, ITAE de los modelos disponibles ===
% Referencia: salida_filt_plot evaluada en la malla uniforme t_u

% 0) Asegurar malla uniforme y referencia yu
if ~(exist('t_u','var'))
    % crear malla uniforme desde tus vectores ya graficados
    tt = tiempo_plot(:); yy = salida_filt_plot(:);
    m  = isfinite(tt) & isfinite(yy);
    tt = tt(m); yy = yy(m);
    [tt, ia] = unique(tt,'stable'); yy = yy(ia);
    dt = median(diff(tt)); if ~isfinite(dt) || dt<=0, dt = mean(diff(tt)); end
    t_u = (tt(1):dt:tt(end)).';
end
yu = interp1(tiempo_plot(:), salida_filt_plot(:), t_u, 'linear','extrap');  % referencia
t_rel = t_u - t_u(1);  % tiempo relativo para ITAE

% 1) Helper para integrar por trapecios
compute_metrics = @(yhat) struct( ...
    'IAE',  trapz(t_u, abs(yu - yhat)), ...
    'ISE',  trapz(t_u, (yu - yhat).^2), ...
    'ITAE', trapz(t_u, t_rel .* abs(yu - yhat)) );

% 2) Recolectar salidas de cada modelo (interpoladas a t_u si hace falta)
models = struct();

if exist('y_model','var') && ~isempty(y_model)
    models.ALFARO_123C = y_model(:);
end
if exist('y_ho','var') && ~isempty(y_ho)
    models.HO = y_ho(:);
end
if exist('y_p1d','var') && ~isempty(y_p1d)
    models.P1D_TOOLBOX = y_p1d(:);
end
if exist('y_LS','var') && ~isempty(y_LS)
    % y_LS se generó sobre 't' (no necesariamente t_u); interpólalo
    yls = interp1(tiempo_plot(:), y_LS(:), t_u, 'linear','extrap');
    models.MINIMOS_CUADRADOS = yls(:);
end

% 3) Calcular métricas y mostrar tabla
fprintf('\n=== Índices de desempeño (ref: salida filtrada) ===\n');
fprintf('%-22s  %12s  %12s  %12s\n', 'Modelo', 'IAE', 'ISE', 'ITAE');
fprintf('%s\n', repmat('-',1,62));

model_names = fieldnames(models);
for k = 1:numel(model_names)
    name = model_names{k};
    yhat = models.(name);
    % asegurar longitudes coincidentes
    if numel(yhat) ~= numel(t_u)
        yhat = interp1(tiempo_plot(:), yhat, t_u, 'linear','extrap');
    end
    met = compute_metrics(yhat);
    fprintf('%-22s  %12.6g  %12.6g  %12.6g\n', name, met.IAE, met.ISE, met.ITAE);
end
fprintf('%s\n', repmat('-',1,62));

%% === Imprimir función de transferencia del modelo P1D (Toolbox) ===
if exist('P1D','var')
    fprintf('\n=== Función de Transferencia del Modelo P1D (Toolbox) ===\n');

    % Si es idproc, convertirlo a tf estándar
    if isa(P1D,'idproc')
        G_P1D_tf = tf(idtf(P1D));
        Kp = P1D.Kp;
        Tp = P1D.Tp1;
        Td = P1D.Td;
    elseif isa(P1D,'idtf')
        G_P1D_tf = tf(P1D);
        Kp = dcgain(G_P1D_tf);
        % obtener parámetros aproximados si es posible
        Tp = 1 / abs(pole(G_P1D_tf));
        Td = G_P1D_tf.InputDelay;
    elseif isa(P1D,'tf')
        G_P1D_tf = P1D;
        Kp = dcgain(G_P1D_tf);
        Tp = 1 / abs(pole(G_P1D_tf));
        Td = G_P1D_tf.InputDelay;
    else
        error('El modelo P1D no es de tipo idproc, idtf o tf.');
    end

    % Imprimir parámetros principales
    fprintf('Ganancia (Kp): %.6g\n', Kp);
    fprintf('Constante de tiempo (Tp): %.6g s\n', Tp);
    fprintf('Tiempo muerto (Td): %.6g s\n', Td);

    % Imprimir la función de transferencia en formato simbólico
    fprintf('\nG_P1D(s) = %.6g / (%.6g·s + 1) · e^{−%.6g·s}\n', Kp, Tp, Td);
    disp(' ');
    disp(G_P1D_tf);
else
    warning('No se encontró el modelo P1D en el workspace.');
end


%% === GRAFICOS FINALES SIMPLIFICADOS CON FIT EN LEYENDA ===

%% FIGURA 1: Referencia + salida filtrada + modelos de identificación
figure;
hold on;

% Referencia (entrada normalizada)


% Señal filtrada (salida real usada como "verdad")
plot(tiempo_plot, salida_filt_plot, 'r-', 'LineWidth', 2, ...
    'DisplayName', 'Salida filtrada');

% Modelo 123c (Alfaro)
if exist('y_model','var') && exist('t_u','var')
    if exist('fit_123c','var')
        lbl_123c = sprintf('Modelo 123c (Fit %.2f%%)', fit_123c);
    else
        lbl_123c = 'Modelo 123c';
    end
    plot(t_u, y_model, 'b-', 'LineWidth', 1.3, ...
        'DisplayName', lbl_123c);
end

% Modelo Ho et al.
if exist('y_ho','var') && exist('t_u','var')
    if exist('fit_ho','var')
        lbl_ho = sprintf('Ho et al. (Fit %.2f%%)', fit_ho);
    else
        lbl_ho = 'Ho et al.';
    end
    plot(t_u, y_ho, 'm-.', 'LineWidth', 1.3, ...
        'DisplayName', lbl_ho);
end

% Modelo P1D (Toolbox)
if exist('y_p1d','var') && exist('t_u','var')
    if exist('fit_txt','var')
        % fit_txt ya viene como string con 2 decimales: '95.23'
        lbl_p1d = sprintf('Modelo P1D (Fit %s%%)', fit_txt);
    else
        lbl_p1d = 'Modelo P1D';
    end
    plot(t_u, y_p1d, 'g-', 'LineWidth', 1.3, ...
        'DisplayName', lbl_p1d);
end

% Modelo por mínimos cuadrados (ARX 1,1,1)
if exist('y_LS','var') && exist('t','var')
    if exist('fit_LS','var')
        lbl_ls = sprintf('Mínimos cuadrados (Fit %.2f%%)', fit_LS);
    else
        lbl_ls = 'Mínimos cuadrados';
    end
    plot(t, y_LS, 'c-', 'LineWidth', 1.3, ...
        'DisplayName', lbl_ls);
end

% Marcar t = 0 pero SIN entrar en la leyenda
xline(0, 'k:', 'HandleVisibility', 'off');
yline(0, 'k-', 'HandleVisibility', 'off');

grid on;
xlabel('Tiempo [s]');
ylabel('Salida / Referencia (normalizadas)');
title('Referencia, salida filtrada y modelos de identificación');
legend('Location', 'southeast');
hold off;



% figura             
%   - SOLO 38–43 s
%   - Muestra señal de control 90→120 y salida real


mask_f1 = (tiempo >= 38) & (tiempo <= 43);

figure; clf;

yyaxis left;
plot(tiempo(mask_f1), senal_control(mask_f1), 'b--', 'LineWidth', 1.8);
ylabel('Control (PWM)');
ylim([85 125]);     % Se observa 90→120 claramente

yyaxis right;
plot(tiempo(mask_f1), salida_velocidad(mask_f1), 'r-', 'LineWidth', 2);
ylabel('Salida (rpm)');

% Colores de ejes en negro (sin colores laterales)
ax = gca;
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = [0 0 0];

xlabel('Tiempo (s)');
xlim([38 43]);
grid on;
title('Datos reales del motor en función del tiempo(s)');
legend('Control','Salida','Location','southeast');

%% === FIGURA 2: Alinear referencia con 40 rpm del eje derecho ===

entrada_ref_30 = 30 * (tiempo_plot >= 0);

% Normalizar salida filtrada
yss_filt = mean(salida_filt_plot(end-20:end));
if isnan(yss_filt), yss_filt = salida_filt_plot(end); end
escala = 30 / yss_filt;
salida_norm_30 = salida_filt_plot * escala;

figure; clf;

%% === EJE IZQUIERDO ===
yyaxis left
plot(tiempo_plot, entrada_ref_30, 'r--', 'LineWidth', 1.8, 'DisplayName','Referencia');
hold on
plot(tiempo_plot, salida_norm_30, 'b-', 'LineWidth', 1.8, 'DisplayName','Salida filtrada (normalizada)');
ylabel('Control');

%% === EJE DERECHO ===
yyaxis right
plot(tiempo(mask_f1), salida_velocidad(mask_f1), 'k-', 'LineWidth', 2);
ylabel('Salida y(s)');
ylim([-10 50]); % límite derecho

%% === ALINEACIÓN DE 30 PWM con 40 rpm ===
% Alturas relativas
y_right_min = -10;
y_right_max = 50;
target_rpm = 40;

fraction = (target_rpm - y_right_min) / (y_right_max - y_right_min);
left_min = 0;
left_max = entrada_ref_30( end );  % normalmente 30

% nuevo límite superior del eje izquierdo
new_left_max = left_max / fraction;
yyaxis left
ylim([left_min new_left_max]);

%% Estética
ax = gca;
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = [0 0 0];
xline(0,'k:');
grid on;
xlabel('Tiempo (s)');
title('Datos filtrados en función del tiempo(s) sin offset ');
legend('Location','southeast');


%% === Cortar la gráfica EXACTAMENTE donde termina la señal filtrada ===

% Último índice válido de la señal filtrada
idx_final = find(~isnan(salida_norm_30), 1, 'last');

% Tiempo correspondiente al último dato
t_final = tiempo_plot(idx_final);

% Limitar el eje X exactamente hasta ahí
xlim([tiempo_plot(1), t_final]);

hold off;
