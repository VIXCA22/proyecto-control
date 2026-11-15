%% === Comparación general: referencia, perturbación y salidas (líneas sólidas) ===
referencia   = out.referencia;
perturbacion = out.perturbacion;

x0 = out.salida;   % Chien
x1 = out.salida1;  % Síntesis analítica
x2 = out.salida2;  % Fertik y Sharpe
x3 = out.salida3;  % LGR discreto

% === SEÑALES DE CONTROL ===
c0 = out.control;    % Chien
c1 = out.control1;   % Síntesis analítica
c2 = out.control2;   % Fertik y Sharpe
c3 = out.control3;   % LGR discreto

u0 = c0.Data;
u1 = c1.Data;
u2 = c2.Data;
u3 = c3.Data;

t = referencia.Time;
r = referencia.Data;
p = perturbacion.Data;

y0 = x0.Data;
y1 = x1.Data;
y2 = x2.Data;
y3 = x3.Data;

%% === Figura comparativa ===
figure; hold on; grid on;
color_ref  = [0.00 0.50 0.00]; 
color_pert = [0.60 0.00 0.60]; 
color_y0   = [0.00 0.45 0.74];  
color_y1   = [0.85 0.33 0.10];  
color_y2   = [0.47 0.67 0.19];  
color_y3   = [1.0 0.85 0.0]; 

plot(t, r,  'Color', color_ref,  'LineWidth', 2, 'DisplayName','Referencia (escalón)');
plot(t, p,  'Color', color_pert, 'LineWidth', 2, 'DisplayName','Perturbación');
plot(t, y0, 'Color', color_y0,   'LineWidth', 2, 'DisplayName','Método Chien(1988)');
plot(t, y1, 'Color', color_y1,   'LineWidth', 2, 'DisplayName','Sintonización analítica');
plot(t, y2, 'Color', color_y2,   'LineWidth', 2, 'DisplayName','Fertik and Sharpe(1979)');
plot(t, y3, 'Color', color_y3,   'LineWidth', 2, 'DisplayName','LGR discreto');
xlabel('Tiempo (s)'); ylabel('Amplitud (V)');
title('Comparación de señales del sistema continuo');
legend('Location','best'); hold off;

%% === Detección robusta del instante del escalón en la referencia (t0) ===
dt   = median(diff(t));
drdt = [0; diff(r)./max(diff(t),eps)];
[~, kmax] = max(abs(drdt));

r0   = median(r(1:max(5,round(0.1*numel(r)))),'omitnan');
rss  = median(r(end-max(5,round(0.1*numel(r)))+1:end),'omitnan');
rmid = r0 + 0.5*(rss - r0);

k2   = find( (r(1:end-1)-rmid).*(r(2:end)-rmid) <= 0, 1, 'first');
if isempty(k2), k2 = kmax; end

t1 = t(k2); 
t2 = t(min(k2+1,numel(t)));
r1 = r(k2); 
r2 = r(min(k2+1,numel(r)));

alpha = (rmid - r1) / max(r2 - r1, eps);
t0 = t1 + alpha*(t2 - t1);

%% === Métricas generales por salida: Mpn(global)%, Ts5(global), Emax ===
metrics = @(y) local_metrics(y, t, t0);

M0 = metrics(y0);   % Chien
M1 = metrics(y1);   % Síntesis
M2 = metrics(y2);   % Fertik
M3 = metrics(y3);   % LGR

fprintf('\n=== Métricas globales (normalizadas a Δy = y_{ss} - y_0) ===\n');
fprintf('%-12s | %-10s | %-10s | %-10s | %-10s | %-10s\n', ...
        'Sistema','y0','yss','Δy','Mpn(%)','Ts5% [s]');
fprintf('-----------------------------------------------------------------------\n');
fprintf('%-12s | %10.4f | %10.4f | %10.4f | %10.3f | %10.4f\n', ...
        'Salida base', M0.y0, M0.yss, M0.dy, M0.Mpn, M0.Ts5);
fprintf('%-12s | %10.4f | %10.4f | %10.4f | %10.3f | %10.4f\n', ...
        'Salida 1',    M1.y0, M1.yss, M1.dy, M1.Mpn, M1.Ts5);
fprintf('%-12s | %10.4f | %10.4f | %10.4f | %10.3f | %10.4f\n', ...
        'Salida 2',    M2.y0, M2.yss, M2.dy, M2.Mpn, M2.Ts5);
fprintf('%-12s | %10.4f | %10.4f | %10.4f | %10.3f | %10.4f\n', ...
        'Salida 3',    M3.y0, M3.yss, M3.dy, M3.Mpn, M3.Ts5);
fprintf('-----------------------------------------------------------------------\n');
fprintf('%-12s | %s\n','Emax (base)', sprintf('%.4f', M0.Emax));
fprintf('%-12s | %s\n','Emax (sal1)', sprintf('%.4f', M1.Emax));
fprintf('%-12s | %s\n','Emax (sal2)', sprintf('%.4f', M2.Emax));
fprintf('%-12s | %s\n','Emax (sal3)', sprintf('%.4f', M3.Emax));

%% === ÍNDICES SERVOCONTROL (0–3.3 s) ===
idx = (t >= 0) & (t <= 3.3);

tt     = t(idx);
r_seg  = r(idx);
y0_seg = y0(idx);
y1_seg = y1(idx);
y2_seg = y2(idx);
y3_seg = y3(idx);

u0_seg = u0(idx);
u1_seg = u1(idx);
u2_seg = u2(idx);
u3_seg = u3(idx);

tol = 1e-3;

fprintf('\n============================================\n');
fprintf('   ÍNDICES SERVOCONTROL (0–3.3 s)\n');
fprintf('   IAEr, TVur, Umax, ta5\n');
fprintf('============================================\n\n');

fprintf('%-22s | %10s | %10s | %10s | %10s\n', ...
        'Controlador','IAEr','TVur','Umax','ta5 (s)');
fprintf('--------------------------------------------------------------------------\n');

% ----- CHIEN (1988) -----
e0    = r_seg - y0_seg;
IAEr0 = trapz(tt, abs(e0));
TVur0 = sum(abs(diff(u0_seg)));
Umax0 = max(abs(max(u0_seg)), abs(min(u0_seg)));

y0_ini = y0_seg(1);   
y0_fin = y0_seg(end);
dy0    = y0_fin - y0_ini;
y0_95  = y0_ini + 0.95*dy0;
y0_105 = y0_ini + 1.05*dy0;

ind0_95  = find(abs(y0_seg - y0_95)  < tol, 1, 'last');
ind0_105 = find(abs(y0_seg - y0_105) < tol, 1, 'last');
ind0_ta5 = max([ind0_95, ind0_105]);
ta5_0    = tt(ind0_ta5) - tt(1);

fprintf('%-22s | %10.4f | %10.4f | %10.4f | %10.4f\n', ...
        'Chien (1988)', IAEr0, TVur0, Umax0, ta5_0);

% ----- SÍNTESIS ANALÍTICA -----
e1    = r_seg - y1_seg;
IAEr1 = trapz(tt, abs(e1));
TVur1 = sum(abs(diff(u1_seg)));
Umax1 = max(abs(max(u1_seg)), abs(min(u1_seg)));

y1_ini = y1_seg(1);   
y1_fin = y1_seg(end);
dy1    = y1_fin - y1_ini;
y1_95  = y1_ini + 0.95*dy1;
y1_105 = y1_ini + 1.05*dy1;

ind1_95  = find(abs(y1_seg - y1_95)  < tol, 1, 'last');
ind1_105 = find(abs(y1_seg - y1_105) < tol, 1, 'last');
ind1_ta5 = max([ind1_95, ind1_105]);
ta5_1    = tt(ind1_ta5) - tt(1);

fprintf('%-22s | %10.4f | %10.4f | %10.4f | %10.4f\n', ...
        'S. Analítica', IAEr1, TVur1, Umax1, ta5_1);

% ----- FERTIK Y SHARPE -----
e2    = r_seg - y2_seg;
IAEr2 = trapz(tt, abs(e2));
TVur2 = sum(abs(diff(u2_seg)));
Umax2 = max(abs(max(u2_seg)), abs(min(u2_seg)));

y2_ini = y2_seg(1);   
y2_fin = y2_seg(end);
dy2    = y2_fin - y2_ini;
y2_95  = y2_ini + 0.95*dy2;
y2_105 = y2_ini + 1.05*dy2;

ind2_95  = find(abs(y2_seg - y2_95)  < tol, 1, 'last');
ind2_105 = find(abs(y2_seg - y2_105) < tol, 1, 'last');
ind2_ta5 = max([ind2_95, ind2_105]);
ta5_2    = tt(ind2_ta5) - tt(1);

fprintf('%-22s | %10.4f | %10.4f | %10.4f | %10.4f\n', ...
        'Fertik-Sharpe', IAEr2, TVur2, Umax2, ta5_2);

% ----- LGR DISCRETO -----
e3    = r_seg - y3_seg;
IAEr3 = trapz(tt, abs(e3));
TVur3 = sum(abs(diff(u3_seg)));
Umax3 = max(abs(max(u3_seg)), abs(min(u3_seg)));

y3_ini = y3_seg(1);   
y3_fin = y3_seg(end);
dy3    = y3_fin - y3_ini;
y3_95  = y3_ini + 0.95*dy3;
y3_105 = y3_ini + 1.05*dy3;

ind3_95  = find(abs(y3_seg - y3_95)  < tol, 1, 'last');
ind3_105 = find(abs(y3_seg - y3_105) < tol, 1, 'last');
ind3_ta5 = max([ind3_95, ind3_105]);
ta5_3    = tt(ind3_ta5) - tt(1);

fprintf('%-22s | %10.4f | %10.4f | %10.4f | %10.4f\n', ...
        'LGR Discreto', IAEr3, TVur3, Umax3, ta5_3);

fprintf('--------------------------------------------------------------------------\n');

%% === SOBREPASO SOLO EN LA VENTANA DE SERVO (0–3.3 s) ===
calc_Mpn_servo = @(tt, yseg) local_Mpn_servo(tt, yseg);

Mpn0_servo = calc_Mpn_servo(tt, y0_seg);
Mpn1_servo = calc_Mpn_servo(tt, y1_seg);
Mpn2_servo = calc_Mpn_servo(tt, y2_seg);
Mpn3_servo = calc_Mpn_servo(tt, y3_seg);

%% === CUADRO II: Índices de bondad para servocontrol (0–3.3 s) ===
fprintf('\n============================================\n');
fprintf('           CUADRO II – SERVOCONTROL\n');
fprintf('  IAEr – ISEr – Mpn%% – Tuwr%%  (sin Ta2%%)\n');
fprintf('============================================\n\n');

fprintf('%-22s | %10s | %10s | %10s | %10s\n', ...
        'Controlador','IAEr','ISEr','Mpn%','Tuwr%');
fprintf('--------------------------------------------------------------------------\n');

% ISEr NORMALIZADO (error dividido por Δy global de cada salida)
E0_norm = e0 ./ max(abs(M0.dy), eps);
ISEr0   = trapz(tt, E0_norm.^2);
Tuwr0   = sum(abs(diff(u0_seg)));

E1_norm = e1 ./ max(abs(M1.dy), eps);
ISEr1   = trapz(tt, E1_norm.^2);
Tuwr1   = sum(abs(diff(u1_seg)));

E2_norm = e2 ./ max(abs(M2.dy), eps);
ISEr2   = trapz(tt, E2_norm.^2);
Tuwr2   = sum(abs(diff(u2_seg)));

E3_norm = e3 ./ max(abs(M3.dy), eps);
ISEr3   = trapz(tt, E3_norm.^2);
Tuwr3   = sum(abs(diff(u3_seg)));

fprintf('%-22s | %10.4f | %10.4f | %10.3f | %10.3f\n', ...
        'Chien (1988)',       IAEr0, ISEr0, Mpn0_servo, Tuwr0);
fprintf('%-22s | %10.4f | %10.4f | %10.3f | %10.3f\n', ...
        'Síntesis Analítica', IAEr1, ISEr1, Mpn1_servo, Tuwr1);
fprintf('%-22s | %10.4f | %10.4f | %10.3f | %10.3f\n', ...
        'Fertik y Sharpe',    IAEr2, ISEr2, Mpn2_servo, Tuwr2);
fprintf('%-22s | %10.4f | %10.4f | %10.3f | %10.3f\n', ...
        'LGR Discreto',       IAEr3, ISEr3, Mpn3_servo, Tuwr3);

fprintf('--------------------------------------------------------------------------\n');

%% === CUADRO III: Índices de bondad para CONTROL REGULATORIO (9–15 s) ===
fprintf('\n============================================\n');
fprintf('           CUADRO III – REGULATORIO\n');
fprintf('     IAEd – ISEd – Tuwd%% – Emax  (9–15 s)\n');
fprintf('============================================\n\n');

fprintf('%-22s | %10s | %10s | %10s | %10s\n', ...
        'Controlador','IAEd','ISEd','Tuwd%','Emax');
fprintf('--------------------------------------------------------------------------\n');

idxr = (t >= 9) & (t <= 15);

tt_r  = t(idxr);
y0_r  = y0(idxr);
y1_r  = y1(idxr);
y2_r  = y2(idxr);
y3_r  = y3(idxr);

u0_r = u0(idxr);
u1_r = u1(idxr);
u2_r = u2(idxr);
u3_r = u3(idxr);

% Nivel estacionario previo a la perturbación (8.5–9 s)
idx_ss = (t >= 8.5 & t <= 9);
y0_ss = mean(y0(idx_ss));
y1_ss = mean(y1(idx_ss));
y2_ss = mean(y2(idx_ss));
y3_ss = mean(y3(idx_ss));

% Errores regulatorios (respecto a y_ss previo)
e0d = y0_r - y0_ss;
e1d = y1_r - y1_ss;
e2d = y2_r - y2_ss;
e3d = y3_r - y3_ss;

IAEd0 = trapz(tt_r, abs(e0d));
ISEd0 = trapz(tt_r, e0d.^2);
Tuwd0 = sum(abs(diff(u0_r)));
Emax0 = max(abs(e0d));

fprintf('%-22s | %10.4f | %10.4f | %10.4f | %10.4f\n', ...
        'Chien (1988)', IAEd0, ISEd0, Tuwd0, Emax0);

IAEd1 = trapz(tt_r, abs(e1d));
ISEd1 = trapz(tt_r, e1d.^2);
Tuwd1 = sum(abs(diff(u1_r)));
Emax1 = max(abs(e1d));

fprintf('%-22s | %10.4f | %10.4f | %10.4f | %10.4f\n', ...
        'Síntesis Analítica', IAEd1, ISEd1, Tuwd1, Emax1);

IAEd2 = trapz(tt_r, abs(e2d));
ISEd2 = trapz(tt_r, e2d.^2);
Tuwd2 = sum(abs(diff(u2_r)));
Emax2 = max(abs(e2d));

fprintf('%-22s | %10.4f | %10.4f | %10.4f | %10.4f\n', ...
        'Fertik y Sharpe', IAEd2, ISEd2, Tuwd2, Emax2);

IAEd3 = trapz(tt_r, abs(e3d));
ISEd3 = trapz(tt_r, e3d.^2);
Tuwd3 = sum(abs(diff(u3_r)));
Emax3 = max(abs(e3d));

fprintf('%-22s | %10.4f | %10.4f | %10.4f | %10.4f\n', ...
        'LGR Discreto', IAEd3, ISEd3, Tuwd3, Emax3);

fprintf('--------------------------------------------------------------------------\n');

%% === RESUMEN Ts5 Y SOBREPASO MÁXIMO ===
fprintf('\n============================================\n');
fprintf('   RESUMEN Ts5 Y SOBREPASO MÁXIMO\n');
fprintf('   (Servo y Regulatorio)\n');
fprintf('============================================\n\n');

Ts5_servo_0 = M0.Ts5;
Ts5_servo_1 = M1.Ts5;
Ts5_servo_2 = M2.Ts5;
Ts5_servo_3 = M3.Ts5;

% Mpn de servo ya corregido (solo 0–3.3 s)
% Mpn0_servo..Mpn3_servo arriba

% Ts5 de regulatorio: error e_d(t) respecto a 0 (rechazo completo)
calc_Ts5_reg = @(tt, e) local_Ts5_error(tt, e);

Ts5_reg_0 = calc_Ts5_reg(tt_r, e0d);
Ts5_reg_1 = calc_Ts5_reg(tt_r, e1d);
Ts5_reg_2 = calc_Ts5_reg(tt_r, e2d);
Ts5_reg_3 = calc_Ts5_reg(tt_r, e3d);

fprintf('%-22s | %12s | %12s | %12s\n', ...
        'Controlador','Mpn servo (%)','Ts5 servo (s)','Ts5 reg (s)');
fprintf('---------------------------------------------------------------------\n');
fprintf('%-22s | %12.3f | %12.4f | %12.4f\n', ...
        'Chien (1988)',       Mpn0_servo, Ts5_servo_0, Ts5_reg_0);
fprintf('%-22s | %12.3f | %12.4f | %12.4f\n', ...
        'Síntesis Analítica', Mpn1_servo, Ts5_servo_1, Ts5_reg_1);
fprintf('%-22s | %12.3f | %12.4f | %12.4f\n', ...
        'Fertik y Sharpe',    Mpn2_servo, Ts5_servo_2, Ts5_reg_2);
fprintf('%-22s | %12.3f | %12.4f | %12.4f\n', ...
        'LGR Discreto',       Mpn3_servo, Ts5_servo_3, Ts5_reg_3);
fprintf('---------------------------------------------------------------------\n');
fprintf('Nota: en regulatorio se reporta Ts5 del error e_d(t); el sobrepaso\n');
fprintf('no se define de forma clásica porque la salida vuelve al mismo valor.\n');

%% ===== Subfunción: métricas globales (servo + regulatorio mezclado) =====
function R = local_metrics(y, t, t0)
    pre  = t < t0;
    if ~any(pre)
        kPre = max(5, round(0.05*numel(t)));
        pre(1:kPre) = true;
    end

    tail = false(size(t));
    kTail = max(5, round(0.1*numel(t)));
    tail(end-kTail+1:end) = true;

    y0  = mean(y(pre),'omitnan');
    yss = mean(y(tail),'omitnan');
    dy  = yss - y0;

    idx0 = find(t>=t0,1,'first');
    tt = t(idx0:end);
    yy = y(idx0:end);

    sgn = sign(dy); if sgn==0, sgn = 1; end
    Mpn = 100 * max(((yy - yss).*sgn) / max(abs(dy), eps));
    if ~isfinite(Mpn) || Mpn < 0, Mpn = 0; end

    band = 0.05 * abs(dy);
    Ts = NaN;

    inside = abs(yy - yss) <= band;
    k1 = find(inside,1,'first');
    if ~isempty(k1)
        for k = k1:numel(yy)
            if all(abs(yy(k:end) - yss) <= band)
                Ts = tt(k) - t0;
                break;
            end
        end
    end

    Emax = max(abs(yy - yss));

    R = struct('y0',y0,'yss',yss,'dy',dy,'Mpn',Mpn,'Ts5',Ts,'Emax',Emax);
end

%% ===== Subfunción: sobrepaso SOLO en ventana de servo =====
function Mpn = local_Mpn_servo(t, y)
    t = t(:);
    y = y(:);

    if numel(y) < 5
        Mpn = 0;
        return;
    end

    k_ini = max(3, round(0.1*numel(y)));
    y0 = mean(y(1:k_ini), 'omitnan');

    k_fin = max(3, round(0.1*numel(y)));
    yss = mean(y(end-k_fin+1:end), 'omitnan');

    dy = yss - y0;
    if abs(dy) < eps
        Mpn = 0;
        return;
    end

    sgn = sign(dy);
    yy_norm = (y - yss) * sgn / abs(dy);

    Mpn = 100 * max(yy_norm);
    if ~isfinite(Mpn) || Mpn < 0
        Mpn = 0;
    end
end

%% ===== Subfunción: Ts5 del error regulatorio =====
function Ts5 = local_Ts5_error(t, e)
    t  = t(:);
    e  = e(:);

    Apeak = max(abs(e));
    if Apeak == 0
        Ts5 = 0;
        return;
    end

    band = 0.05 * Apeak;
    inside = abs(e) <= band;
    Ts5 = NaN;

    k1 = find(inside, 1, 'first');
    if ~isempty(k1)
        for k = k1:numel(e)
            if all(inside(k:end))
                Ts5 = t(k) - t(1);
                break;
            end
        end
    end
end
