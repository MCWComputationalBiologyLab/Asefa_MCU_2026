close all; clear; clc;

ModelFolder = 'M:\Dash lab\Yisak_Asefa\MCU_Manuscript_Supplement\Simulation\Cation_Fluxes';
addpath(ModelFolder);

%% ----------------- Globals used by the ODE -----------------
global t_Ca_add t_Na_add g_Ca_amplitude g_pulse_width_Ca
global USE_G_UNITY g0 g1 g12
global USE_INACTIVATION

% ---- Toggle P_open_x here ----
% USE_INACTIVATION = 0  -->  P_open_x = 1          (inactivation OFF)
% USE_INACTIVATION = 1  -->  P_open_x = S_O + S_R  (inactivation ON)
USE_INACTIVATION = 1;   % <--- change this line to switch behaviour
% ------------------------------

% Black lines (USE_G_UNITY = 0): MICU gating ON
%   --> always "With Biphasic [Ca2+]_m Regulation of MCU Kinetics"
% Red lines  (USE_G_UNITY = 1): gate forced open
%   --> always "No [Ca2+]_m Regulation of MCU Kinetics"
labelBlack = 'With Biphasic [Ca^{2+}]_m Regulation of MCU Kinetics';
labelRed   = 'No [Ca^{2+}]_m Regulation of MCU Kinetics';

% Common MICU gating weights (must match ODE’s meaning of x0/x1/x12)
g0  = 0.0;
g1  = 0.5;
g12 = 0.95;

% Common pulse timing (as in both scripts)
t_Ca_add         = 700;   % time we add Ca2+ (s)
t_Na_add         = 360;   % (compatibility only; Na pulse not used)
g_pulse_width_Ca = 4;     % Ca pulse width (s)

% Common ODE solver options (as in both scripts)
options = odeset('RelTol',1e-12,'AbsTol',1e-10, ...
                 'MaxStep',0.2, 'NonNegative',1:6, 'Refine',1);

%% =======================================================================
%  RUN 1: Fig6_Only_N7 settings
% =======================================================================
% Initial conditions (N7)
% y = [K_e,   K_x,   Na_e,  Na_x,   Ca_e,   Ca_x,   x0,  x1,  S_O, S_I, S_R]
IC_N7 = [100e-3, 150e-3,  2e-3,   0,   20e-9, 100e-9,  1,   0,   1,   0,   0];

t_end_N7 = 1100;                           % total time (s)
Ca_amplitudes_N7 = [0.3, 0.7, 1.0, 1.5]*1e-6;

% Figure:
figure(1); clf; set(gcf,'Position',[20,120,1500,570],'Color','w'); hold on

axBlank = nexttile; axis(axBlank,'off');

axPopen = nexttile; hold(axPopen,'on');
ylabel(axPopen,'P_{Open} (unitless)'); xlabel(axPopen,'Time (s)');
xlim(axPopen,[650 1000]); xticks(axPopen,700:100:1100); set(axPopen,'FontSize',16);
ylim(axPopen,[0 0.3]); yticks(axPopen,0:0.05:0.16);

axCx = nexttile; hold(axCx,'on');
ylabel(axCx,'[Ca^{2+}]_m (\muM)'); xlabel(axCx,'Time (s)');
xlim(axCx,[650 1000]); xticks(axCx,700:100:1100); set(axCx,'FontSize',16);
ylim(axCx,[0.08 .21]); yticks(axCx,0.0:0.04:0.25);

axCe = nexttile; hold(axCe,'on');
ylabel(axCe,'[Ca^{2+}]_c (nM)'); xlabel(axCe,'Time (s)');
xlim(axCe,[650 1000]); xticks(axCe,700:100:1100); set(axCe,'FontSize',16);
ylim(axCe,[0 2000]); yticks(axCe,0:500:2000);   % (fix: was xticks in your text)

axJmcu = nexttile; hold(axJmcu,'on');
ylabel(axJmcu,'J_{MCU} (nmol mg^{-1} s^{-1})'); xlabel(axJmcu,'Time (s)');
xlim(axJmcu,[650 1100]); set(axJmcu,'FontSize',16);

axBPCa2 = nexttile; hold(axBPCa2,'on');
ylabel(axBPCa2,'[BPCa]_m (\muM)'); xlabel(axBPCa2,'Time (s)');
xlim(axBPCa2,[650 1000]); xticks(axBPCa2,700:100:1100); set(axBPCa2,'FontSize',16);
ylim(axBPCa2,[700 2100]); yticks(axBPCa2,00:500:2000);

% ---- Sweep: N7 ----
for i = 1:numel(Ca_amplitudes_N7)
    g_Ca_amplitude = Ca_amplitudes_N7(i);

    % (A) MICU gating enabled (black)
    USE_G_UNITY = 0;
    [tM, yM] = ode15s(@Cation_Dynamics_ODEs, [0 t_end_N7], IC_N7, options);

    Ca_e_M = yM(:,5);
    Ca_x_M = yM(:,6);
    x0_M   = yM(:,7);
    x1_M   = yM(:,8);
    x12_M  = 1 - x0_M - x1_M;
    G_M    = g0*x0_M + g1*x1_M + g12*x12_M;

    [~, BPCaM, ~] = Ca_Buffering_in(Ca_x_M);
    fluxM = compute_flux_series(tM, yM, 0, g0, g1, g12);

    if i == 1
        plot(axPopen, tM, G_M,            '-', 'Color','k', 'LineWidth',1.1, 'DisplayName', labelBlack);
    else
        plot(axPopen, tM, G_M,            '-', 'Color','k', 'LineWidth',1.1, 'HandleVisibility','off');
    end
    plot(axCx,    tM, Ca_x_M*1e6,     '-', 'Color','k', 'LineWidth',1.3);
    plot(axCe,    tM, Ca_e_M*1e9,     '-', 'Color','k', 'LineWidth',1.1);
    plot(axJmcu,  tM, fluxM.J_MCU*1e3,'-', 'Color','k', 'LineWidth',1.0);
    plot(axBPCa2, tM, BPCaM*1e6,      '-', 'Color','k', 'LineWidth',1.0);

    % (B) Gate forced to 1 (always open, red)
    USE_G_UNITY = 1;
    [tU, yU] = ode15s(@Cation_Dynamics_ODEs, [0 t_end_N7], IC_N7, options);

    Ca_e_U = yU(:,5);
    Ca_x_U = yU(:,6);

    [~, BPCaU, ~] = Ca_Buffering_in(Ca_x_U);
    fluxU = compute_flux_series(tU, yU, 1, g0, g1, g12);

    % Note: P_open not plotted for G=1 (matches your N7 script)
    % Add a dummy invisible line on axPopen for legend on first iteration only
    if i == 1
        plot(axPopen, NaN, NaN, '-', 'Color','r', 'LineWidth',1.1, 'DisplayName', labelRed);
    end
    plot(axCx,    tU, Ca_x_U*1e6,     '-', 'Color','r', 'LineWidth',1.1);
    plot(axCe,    tU, Ca_e_U*1e9,     '-', 'Color','r', 'LineWidth',1.1);
    plot(axJmcu,  tU, fluxU.J_MCU*1e3,'-', 'Color','r', 'LineWidth',1.0);
    plot(axBPCa2, tU, BPCaU*1e6,      '-', 'Color','r', 'LineWidth',1.0);
end

% legend(axPopen, 'show', 'Location', 'northeast', 'FontSize', 11);

%% =======================================================================
%  RUN 2: Fig6_Only settings
% =======================================================================
% Initial conditions (Fig6_Only)
% y = [K_e,   K_x,    Na_e,   Na_x,   Ca_e,   Ca_x,   x0,  x1,  S_O, S_I, S_R]
IC_F6 = [120e-3, 100e-3,  1e-3,  5e-3,  50e-9, 100e-9, 1,   0,   1,   0,   0];

t_end_F6 = 1200;
Ca_amplitudes_F6 = [1.5 2.5 5.0 7.5 10]*1e-6;

% Figure: 
figure(2); clf; set(gcf,'Position',[500,120,1500,570],'Color','w'); hold on

axBlank = nexttile; axis(axBlank,'off');

axPopen = nexttile; hold(axPopen,'on');
ylabel(axPopen,'P_{Open} (unitless)'); xlabel(axPopen,'Time (s)');
xlim(axPopen,[650 1000]); xticks(axPopen,700:100:1100);
ylim(axPopen,[0 1.1]); yticks(axPopen,0:0.2:1); set(axPopen,'FontSize',16);

axCx = nexttile; hold(axCx,'on');
ylabel(axCx,'[Ca^{2+}]_m (\muM)'); xlabel(axCx,'Time (s)');
xlim(axCx,[650 1000]); xticks(axCx,700:100:1100);
ylim(axCx,[0.05 0.55]); yticks(axCx,0:0.1:1.2); set(axCx,'FontSize',16);

axCe = nexttile; hold(axCe,'on');
ylabel(axCe,'[Ca^{2+}]_c (\muM)'); xlabel(axCe,'Time (s)');
xlim(axCe,[650 1000]); xticks(axCe,700:100:1100);
yticks(axCe,0:2:8); set(axCe,'FontSize',16);

axJmcu = nexttile; hold(axJmcu,'on');
ylabel(axJmcu,'J_{MCU} (nmol mg^{-1} s^{-1})'); xlabel(axJmcu,'Time (s)');
xlim(axJmcu,[650 1000]); set(axJmcu,'FontSize',16);

axBPCa2 = nexttile; hold(axBPCa2,'on');
ylabel(axBPCa2,'[BPCa]_m (\muM)'); xlabel(axBPCa2,'Time (s)');
xlim(axBPCa2,[650 1000]); xticks(axBPCa2,700:100:1100); set(axBPCa2,'FontSize',16);
ylim(axBPCa2,[500 5000]); yticks(axBPCa2,500:1500:5000); set(axBPCa2,'FontSize',16);

% ---- Sweep: Fig6_Only ----
for i = 1:numel(Ca_amplitudes_F6)
    g_Ca_amplitude = Ca_amplitudes_F6(i);

    % (A) MICU gating enabled (black)
    USE_G_UNITY = 0;
    [tM, yM] = ode15s(@Cation_Dynamics_ODEs, [0 t_end_F6], IC_F6, options);

    Ca_e_M = yM(:,5);
    Ca_x_M = yM(:,6);
    x0_M   = yM(:,7);
    x1_M   = yM(:,8);
    x12_M  = 1 - x0_M - x1_M;
    G_M    = g0*x0_M + g1*x1_M + g12*x12_M;

    [~, BPCaM, ~] = Ca_Buffering_in(Ca_x_M);
    fluxM = compute_flux_series(tM, yM, 0, g0, g1, g12);

    if i == 1
        plot(axPopen, tM, G_M,            '-', 'Color','k', 'LineWidth',1.1, 'DisplayName', labelBlack);
    else
        plot(axPopen, tM, G_M,            '-', 'Color','k', 'LineWidth',1.1, 'HandleVisibility','off');
    end
    plot(axCx,    tM, Ca_x_M*1e6,     '-', 'Color','k', 'LineWidth',1.3);
    plot(axCe,    tM, Ca_e_M*1e6,     '-', 'Color','k', 'LineWidth',1.1);
    plot(axJmcu,  tM, fluxM.J_MCU*1e3,'-', 'Color','k', 'LineWidth',1.0);
    plot(axBPCa2, tM, BPCaM*1e6,      '-', 'Color','k', 'LineWidth',1.0);

    % (B) Gate forced to 1 (always open, red)
    USE_G_UNITY = 1;
    [tU, yU] = ode15s(@Cation_Dynamics_ODEs, [0 t_end_F6], IC_F6, options);

    Ca_e_U = yU(:,5);
    Ca_x_U = yU(:,6);

    [~, BPCaU, ~] = Ca_Buffering_in(Ca_x_U);
    fluxU = compute_flux_series(tU, yU, 1, g0, g1, g12);

    % --- draw P_open = 1 reference line (red) ---
    if i == 1
        yline(axPopen, 1, '-', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', labelRed);
    else
        yline(axPopen, 1, '-', 'Color', 'r', 'LineWidth', 1.5, 'HandleVisibility','off');
    end

    plot(axCx,    tU, Ca_x_U*1e6,     '-', 'Color','r', 'LineWidth',1.1);
    plot(axCe,    tU, Ca_e_U*1e6,     '-', 'Color','r', 'LineWidth',1.1);
    plot(axJmcu,  tU, fluxU.J_MCU*1e3,'-', 'Color','r', 'LineWidth',1.0);
    plot(axBPCa2, tU, BPCaU*1e6,      '-', 'Color','r', 'LineWidth',1.0);
end

% legend(axPopen, 'show', 'Location', 'northeast', 'FontSize', 11);

%% ========================= Local helper =========================
function flux = compute_flux_series(t, Y, force_open_flag, g0, g1, g12)
% Recompute fluxes from simulated states to match ODE definitions.

% ---- Constants (must mirror ODE) ----
pH_e = 7.2; pH_x = 7.5;
H_e  = 10^-pH_e; H_x = 10^-pH_x;

V_max_NCE = 1.76e-4;
V_max_KHE = 1;
V_max_MKU = 1;
V_max_NHE = 12;
V_max_CHE = 63.96;

Mg_x = 1e-3; %#ok<NASGU>  % kept for mirroring; not used below
Mg_e = 0.5e-3;
Pi_x = 0*2.5e-3; %#ok<NASGU>
Pi_e = 0*2.5e-3; %#ok<NASGU>
dPsi = 180e-3;

x_scale  = 0.05;  % the "x" used in J_MCU in the ODE
G_forced = 0.95;  % matches the ODE's forced-open G

% ---- Unpack states (vectors) ----
K_e  = Y(:,1);  K_x = Y(:,2);
Na_e = Y(:,3);  Na_x = Y(:,4);
Ca_e = Y(:,5);  Ca_x = Y(:,6);
x0   = Y(:,7);  x1   = Y(:,8);
x12  = 1 - x0 - x1;

% ---- Compute G_MICU trace ----
if force_open_flag==0
    G_MICU = g0.*x0 + g1.*x1 + g12.*x12;
else
    G_MICU = G_forced.*ones(size(x0));
end

% ---- Allocate outputs ----
n = numel(t);
J_MKU = zeros(n,1); J_KHE = J_MKU; J_NHE = J_MKU;
J_NCE = J_MKU;      J_CHE = J_MKU; J_MCU = J_MKU;

% ---- Loop (flux functions typically scalar-valued) ----
for k = 1:n
    J_MKU(k) = (0.01) * V_max_MKU * MKU(K_x(k), K_e(k));
    J_KHE(k) = (1)    * V_max_KHE * KHE(K_x(k), K_e(k), H_x, H_e);
    J_NHE(k) = (0.1)  * V_max_NHE * NHE(Na_x(k), Na_e(k), H_x, H_e);
    J_NCE(k) = (0.1)  * V_max_NCE * NCE(Na_x(k), Na_e(k), Ca_x(k), Ca_e(k));
    J_CHE(k) = (0.001)* V_max_CHE * (1/60) * CHE(Ca_x(k), Ca_e(k), H_x, H_e);
    J_MCU(k) = (x_scale * 2 * G_MICU(k)) * MCU(Ca_x(k), Ca_e(k), Mg_e, dPsi);
end

flux.J_MKU = J_MKU; flux.J_KHE = J_KHE; flux.J_NHE = J_NHE;
flux.J_NCE = J_NCE; flux.J_CHE = J_CHE; flux.J_MCU = J_MCU;
end