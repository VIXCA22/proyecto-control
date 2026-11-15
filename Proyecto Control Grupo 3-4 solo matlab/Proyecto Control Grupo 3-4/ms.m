%% Comparación de 3 PI sobre planta de 1er orden: Ms, márgenes y Bode
clear; clc; close all;

%% Parámetros de la planta
K = 1.3247;
T = 0.2963;
s = tf('s');
G = K / (T*s + 1);

%% Controladores PI (forma estándar: Kp * (Ti s + 1)/(Ti s))
% 1) Síntesis analítica
Kp_ana = 0.6176;
Ti_ana = 0.2933;
C_ana  = Kp_ana * (Ti_ana*s + 1)/(Ti_ana*s);

% 2) Método de Chien
Kp_ch  = 0.7549;
Ti_ch  = 0.2963;
C_ch   = Kp_ch * (Ti_ch*s + 1)/(Ti_ch*s);

% 3) Método de Fertik–Sharpe
Kp_fe  = 0.4227;
Ti_fe  = 0.1926;
C_fe   = Kp_fe * (Ti_fe*s + 1)/(Ti_fe*s);

%% Frecuencias de análisis
w = logspace(-3, 3, 4000);   % rad/s

%% Lazos abiertos
L_ana = C_ana * G;
L_ch  = C_ch  * G;
L_fe  = C_fe  * G;

%% Sensibilidades S(s) = 1/(1+L(s))
S_ana = feedback(1, L_ana);
S_ch  = feedback(1, L_ch);
S_fe  = feedback(1, L_fe);

[Smag_ana, ~] = bode(S_ana, w);
[Smag_ch,  ~] = bode(S_ch,  w);
[Smag_fe,  ~] = bode(S_fe,  w);

Smag_ana = squeeze(Smag_ana);
Smag_ch  = squeeze(Smag_ch);
Smag_fe  = squeeze(Smag_fe);

% Ms y frecuencias
[Ms_ana, idx1] = max(Smag_ana);  wMs_ana = w(idx1);
[Ms_ch,  idx2] = max(Smag_ch );  wMs_ch  = w(idx2);
[Ms_fe,  idx3] = max(Smag_fe );  wMs_fe  = w(idx3);

%% Márgenes para cada lazo abierto L = C*G
[Gm_ana, Pm_ana, Wcg_ana, Wcp_ana] = margin(L_ana);
[Gm_ch,  Pm_ch,  Wcg_ch,  Wcp_ch ] = margin(L_ch);
[Gm_fe,  Pm_fe,  Wcg_fe,  Wcp_fe ] = margin(L_fe);

Gm_ana_dB = 20*log10(Gm_ana);
Gm_ch_dB  = 20*log10(Gm_ch);
Gm_fe_dB  = 20*log10(Gm_fe);

%% --- Imprimir resultados en consola ---
fprintf('\n=========== Índice Ms (sensibilidad máxima) ===========\n');
fprintf('Analítica:     Ms = %.4f  en w = %.4f rad/s  (%.2f dB)\n', ...
        Ms_ana, wMs_ana, 20*log10(Ms_ana));
fprintf('Chien:         Ms = %.4f  en w = %.4f rad/s  (%.2f dB)\n', ...
        Ms_ch,  wMs_ch,  20*log10(Ms_ch));
fprintf('Fertik-Sharpe: Ms = %.4f  en w = %.4f rad/s  (%.2f dB)\n', ...
        Ms_fe,  wMs_fe,  20*log10(Ms_fe));

fprintf('\n=========== Márgenes de ganancia y fase ===========\n');
fprintf('Analítica:     GM = %.3f (%.2f dB)  Wcg = %.4f  |  PM = %.2f°  Wcp = %.4f\n', ...
        Gm_ana, Gm_ana_dB, Wcg_ana, Pm_ana, Wcp_ana);
fprintf('Chien:         GM = %.3f (%.2f dB)  Wcg = %.4f  |  PM = %.2f°  Wcp = %.4f\n', ...
        Gm_ch,  Gm_ch_dB,  Wcg_ch,  Pm_ch,  Wcp_ch );
fprintf('Fertik-Sharpe: GM = %.3f (%.2f dB)  Wcg = %.4f  |  PM = %.2f°  Wcp = %.4f\n', ...
        Gm_fe,  Gm_fe_dB,  Wcg_fe,  Pm_fe,  Wcp_fe );

%% ================== BODE LAZO ABIERTO (3 curvas) ==================
[mag_ana, phase_ana] = bode(L_ana, w);
[mag_ch,  phase_ch ] = bode(L_ch,  w);
[mag_fe,  phase_fe ] = bode(L_fe,  w);

mag_ana   = squeeze(mag_ana);
mag_ch    = squeeze(mag_ch);
mag_fe    = squeeze(mag_fe);
phase_ana = squeeze(phase_ana);
phase_ch  = squeeze(phase_ch);
phase_fe  = squeeze(phase_fe);

figure('Name','Bode lazo abierto - 3 PI','NumberTitle','off');
colors = lines(3);

% ---- Magnitud ----
subplot(2,1,1);
semilogx(w, 20*log10(mag_ana), 'Color', colors(1,:), 'LineWidth', 1.8); hold on;
semilogx(w, 20*log10(mag_ch),  'Color', colors(2,:), 'LineWidth', 1.8);
semilogx(w, 20*log10(mag_fe),  'Color', colors(3,:), 'LineWidth', 1.8);
yline(0, '--k');
grid on;
ylabel('Magnitud [dB]');
title('Bode del lazo abierto L(s) - Tres métodos PI');
legend('Analítica','Chien','Fertik-Sharpe','Location','Best');

% Líneas de ωg
if Wcg_ana > 0, xline(Wcg_ana,'--','Color',colors(1,:)); end
if Wcg_ch  > 0, xline(Wcg_ch, '--','Color',colors(2,:)); end
if Wcg_fe  > 0, xline(Wcg_fe, '--','Color',colors(3,:)); end

% Bloque de texto con GM y ωg dentro de la figura (NORMALIZED)
subplot(2,1,1); hold on;
txt1 = sprintf('Analítica:\n  \\omega_g = %.3f\n  GM = %.1f dB', Wcg_ana, Gm_ana_dB);
txt2 = sprintf('Chien:\n  \\omega_g = %.3f\n  GM = %.1f dB', Wcg_ch,  Gm_ch_dB);
txt3 = sprintf('Fertik:\n  \\omega_g = %.3f\n  GM = %.1f dB', Wcg_fe,  Gm_fe_dB);

text(0.72, 0.80, txt1, 'Units','normalized', ...
    'Color', colors(1,:), 'FontSize', 8, 'FontWeight','bold');
text(0.72, 0.55, txt2, 'Units','normalized', ...
    'Color', colors(2,:), 'FontSize', 8, 'FontWeight','bold');
text(0.72, 0.30, txt3, 'Units','normalized', ...
    'Color', colors(3,:), 'FontSize', 8, 'FontWeight','bold');

% ---- Fase ----
subplot(2,1,2);
semilogx(w, phase_ana, 'Color', colors(1,:), 'LineWidth', 1.8); hold on;
semilogx(w, phase_ch,  'Color', colors(2,:), 'LineWidth', 1.8);
semilogx(w, phase_fe,  'Color', colors(3,:), 'LineWidth', 1.8);
yline(-180, '--k');
grid on;
xlabel('\omega [rad/s]');
ylabel('Fase [°]');
legend('Analítica','Chien','Fertik-Sharpe','Location','Best');

% Líneas de ωp
if Wcp_ana > 0, xline(Wcp_ana,'--','Color',colors(1,:)); end
if Wcp_ch  > 0, xline(Wcp_ch, '--','Color',colors(2,:)); end
if Wcp_fe  > 0, xline(Wcp_fe, '--','Color',colors(3,:)); end

% Bloque de texto con PM y ωp dentro de la figura (NORMALIZED)
subplot(2,1,2); hold on;
txt1p = sprintf('Analítica:\n  \\omega_p = %.3f\n  PM = %.1f^\\circ', Wcp_ana, Pm_ana);
txt2p = sprintf('Chien:\n  \\omega_p = %.3f\n  PM = %.1f^\\circ', Wcp_ch,  Pm_ch );
txt3p = sprintf('Fertik:\n  \\omega_p = %.3f\n  PM = %.1f^\\circ', Wcp_fe,  Pm_fe );

text(0.72, 0.80, txt1p, 'Units','normalized', ...
    'Color', colors(1,:), 'FontSize', 8, 'FontWeight','bold');
text(0.72, 0.55, txt2p, 'Units','normalized', ...
    'Color', colors(2,:), 'FontSize', 8, 'FontWeight','bold');
text(0.72, 0.30, txt3p, 'Units','normalized', ...
    'Color', colors(3,:), 'FontSize', 8, 'FontWeight','bold');

%% ================== BODE DE SENSIBILIDAD |S(jw)| (3 curvas + Ms) ==================
figure('Name','Sensibilidad - Ms','NumberTitle','off');
semilogx(w, 20*log10(Smag_ana), 'Color',colors(1,:), 'LineWidth',1.8); hold on;
semilogx(w, 20*log10(Smag_ch ), 'Color',colors(2,:), 'LineWidth',1.8);
semilogx(w, 20*log10(Smag_fe ), 'Color',colors(3,:), 'LineWidth',1.8);
grid on;
xlabel('\omega [rad/s]');
ylabel('|S(j\omega)| [dB]');
title('Sensibilidad S(s) con M_s marcado');
legend('Analítica','Chien','Fertik-Sharpe','Location','Best');

% Marcar puntos Ms
plot(wMs_ana, 20*log10(Ms_ana), 'o', 'Color',colors(1,:), 'MarkerFaceColor',colors(1,:));
plot(wMs_ch,  20*log10(Ms_ch ), 'o', 'Color',colors(2,:), 'MarkerFaceColor',colors(2,:));
plot(wMs_fe,  20*log10(Ms_fe ), 'o', 'Color',colors(3,:), 'MarkerFaceColor',colors(3,:));

text(wMs_ana, 20*log10(Ms_ana)+1, sprintf('Ms=%.2f', Ms_ana), 'Color',colors(1,:), 'FontSize',8);
text(wMs_ch,  20*log10(Ms_ch )+1, sprintf('Ms=%.2f', Ms_ch ), 'Color',colors(2,:), 'FontSize',8);
text(wMs_fe,  20*log10(Ms_fe )+1, sprintf('Ms=%.2f', Ms_fe ), 'Color',colors(3,:), 'FontSize',8);

