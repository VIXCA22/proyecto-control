clear; clc; close all;

%% Planta
K = 1.3247;
T = 0.2963;
s = tf('s');
G = K / (T*s + 1);

%% Controladores PI
% 1) Analítico
Kp_ana = 0.6176;
Ti_ana = 0.2933;
C_ana  = Kp_ana * (Ti_ana*s + 1)/(Ti_ana*s);

% 2) Chien
Kp_ch  = 0.7549;
Ti_ch  = 0.2963;
C_ch   = Kp_ch * (Ti_ch*s + 1)/(Ti_ch*s);

% 3) Fertik–Sharpe
Kp_fe  = 0.4227;
Ti_fe  = 0.1926;
C_fe   = Kp_fe * (Ti_fe*s + 1)/(Ti_fe*s);

controllers = {C_ana, C_ch, C_fe};
names = {'PI Analítico', 'PI Chien', 'PI Fertik–Sharpe'};

%% Rango de frecuencias
w = logspace(-3, 2, 2000);  % 0.001 a 100 rad/s

%% ===== RECORRER LOS 3 CONTROLADORES =====
for k = 1:3
    
    Ck = controllers{k};
    namek = names{k};

    % 1/L(s)
    Lk = Ck * G;
    invL = 1 / Lk;

    % Calcular magnitud, fase y frecuencia
    [mag, phase, wout] = bode(invL, w);
    mag   = squeeze(mag);      % valor absoluto
    phase = squeeze(phase);    % deg
    mag_dB = 20*log10(mag);    % dB

    % Encontrar el punto importante: el máximo de |1/L|
    [peak, idx] = max(mag);
    peak_dB = 20*log10(peak);
    w_peak = wout(idx);

    % ================= FIGURA =================
    figure;

    %% --- MAGNITUD ---
    subplot(2,1,1);
    semilogx(wout, mag_dB, 'LineWidth', 1.5); hold on; grid on;

    % marcar el pico
    plot(w_peak, peak_dB, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    xline(w_peak, ':', 'LineWidth', 1.2);

    text(w_peak, peak_dB, ...
        sprintf('\\omega = %.4f rad/s\n|1/L| = %.4f (%.2f dB)', ...
        w_peak, peak, peak_dB), ...
        'VerticalAlignment','bottom','HorizontalAlignment','left');

    ylabel('Magnitud |1/L(j\\omega)| [dB]');
    title(['Bode de 1/L(s) - ', namek]);

    %% --- FASE ---
    subplot(2,1,2);
    semilogx(wout, phase, 'LineWidth', 1.5); hold on; grid on;
    xline(w_peak, ':', 'LineWidth', 1.2);

    xlabel('Frecuencia (rad/s)');
    ylabel('Fase [^\circ]');
    title('Fase de 1/L(s)');
end
