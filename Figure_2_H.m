
%% -------------------------------------------------------------------------
% MICU gating: time to steady state vs [Ca2+]c using eigenvalues
% For each [Ca2+]c:
%   1) Build the Markov rate matrix Q(Ca).
%   2) Compute its eigenvalues.
%   3) Take the slowest nonzero eigenvalue to get tau_slow.
%   4) Define T_ss = N_tau * tau_slow (e.g. N_tau = 4.6 for ~1%).
% This is exact for linear Markov models and avoids long ODE integrations.
% -------------------------------------------------------------------------

close all; clear; clc;

% [Ca2+] sweep (100 nM -> 100 µM)
C_range = logspace(log10(100e-9), log10(100e-6), 500);  % M

% How close to steady state?  eps = 0.01 -> ~1% from steady state
eps_ss = 0.12;
N_tau  = -log(eps_ss);   % ~4.6

Tss_WT = zeros(size(C_range));
Tss_KO = zeros(size(C_range));

for i = 1:numel(C_range)
    Ca = C_range(i);

    % --------------------------
    % MICU12-WT (3-state) model
    % States: x0 <-> x1 <-> x2
    % --------------------------
    % k1_off = 1.5;                 % MICU1 unbinding (s^-1)
    % k2_off = 2.8;                 % MICU2 unbinding (s^-1)
    % k1_on  = k1_off / (3.5e-6)^2; % MICU1 binding (M^-2 s^-1)
    % k2_on  = k2_off / (5.5e-6)^2; % MICU2 binding (M^-2 s^-1)
    k1_off = 1;                 % MICU1 unbinding (s^-1)
    k2_off = 1;                 % MICU2 unbinding (s^-1)
    k1_on  = 1.1*k1_off / (3.5e-6)^2; % MICU1 binding (M^-2 s^-1)
    k2_on  = 20*k2_off / (9.5e-6)^2; % MICU2 binding (M^-2 s^-1)

    Ca2 = Ca^2;

    % Q_WT: dx/dt = Q_WT * x, with columns summing to zero
    % x = [x0; x1; x2]
    Q_WT = [ ...
        -(k1_on*Ca2),           k1_off,              0          ; ...
         (k1_on*Ca2),  -(k1_off + k2_on*Ca2),      k2_off       ; ...
               0     ,           k2_on*Ca2,       -k2_off       ];

    lam_WT = eig(Q_WT);                 % 3 eigenvalues: 0, λ2, λ3
    lam_WT = real(lam_WT);

    % Remove the zero eigenvalue (steady-state mode)
    lam_nz = lam_WT(abs(lam_WT) > 1e-9);

    % Slowest decay mode: eigenvalue closest to zero (largest real part)
    lam_slow = max(lam_nz);            % negative number
    tau_slow_WT = -1 / lam_slow;       % positive time constant

    Tss_WT(i) = N_tau * tau_slow_WT;

    % --------------------------
    % MICU2-KO (2-state) model
    % States: x0 <-> x1
    % --------------------------
    k1_off_KO = 1.5;
    k1_on_KO  = 1.5 / (1.0e-6)^2;

    % 2x2 rate matrix
    Q_KO = [ ...
        -(k1_on_KO*Ca2),       k1_off_KO; ...
         (k1_on_KO*Ca2),      -k1_off_KO ];

    lam_KO = eig(Q_KO);        % eigenvalues: 0 and -(k_on*Ca^2 + k_off)
    lam_KO = real(lam_KO);
    lam_nz_KO = lam_KO(abs(lam_KO) > 1e-9);
    lam_slow_KO = max(lam_nz_KO);
    tau_slow_KO = -1 / lam_slow_KO;

    Tss_KO(i) = N_tau * tau_slow_KO;
end

% --------------------------
% Plotting: T_ss vs [Ca2+]c
% --------------------------
figure; hold on; box on;

plot(C_range*1e6, Tss_WT, 'k-', 'LineWidth', 2);   % WT
plot(C_range*1e6, Tss_KO, 'b-', 'LineWidth', 2);   % MICU2-KO

xlabel('[Ca^{2+}]_c (\muM)', 'FontSize', 16);
ylabel('T_{ss} (s)', 'FontSize', 16);
legend({'WT (MICU1/2)', 'MICU2-KO'}, 'Location', 'northeast', 'Box', 'off');
xlim([0 20]);
% set(gca, 'XScale', 'log');  % optional if you want log Ca axis
set(gca, 'FontSize', 20, 'LineWidth', 1.5);
