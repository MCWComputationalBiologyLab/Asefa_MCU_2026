close all; clear; clc;

%% Model folder
ModelFolder = 'M:\Dash lab\Yisak_Asefa\MCU_Manuscript_Supplement\Simulation\Cation_Fluxes';
addpath(ModelFolder);
% currentFolder = pwd;
% addpath(fullfile(currentFolder, 'Cation_Fluxes'));

%% Initial conditions
% y = [K_e, K_x, Na_e, Na_x, Ca_e, Ca_x, x0, x1, S_O, S_I, S_R]
IC0 = [100e-3, 150e-3, 2e-3, 0, 100e-9, 100e-9, 1, 0, 1, 0, 0];

%% Globals
global t_Ca_add g_Ca_amplitude g_pulse_width_Ca
global USE_G_UNITY g0 g1 g12
global USE_INACTIVATION

% ---- Toggle P_open_x here ----
% USE_INACTIVATION = 0  -->  P_open_x = 1          (inactivation OFF)
% USE_INACTIVATION = 1  -->  P_open_x = S_O + S_R  (inactivation ON)
USE_INACTIVATION = 1;   % <--- change this line to switch behaviour
% ------------------------------

% Condition label used in figure titles and legend entries
if USE_INACTIVATION
    condLabel = 'With Biphasic [Ca^{2+}]_m Regulation of MCU Kinetics';
else
    condLabel = 'No [Ca^{2+}]_m Regulation of MCU Kinetics';
end

t_end            = 1100;
t_Ca_add         = 700;
g_pulse_width_Ca = 4;

%% WT and KO gating weights
WT_g0  = 0.01;  WT_g1  = 0.45;  WT_g12 = 0.95;
KO_g0  = 1.0;   KO_g1  = 1.0;   KO_g12 = 1.0;

%% Ca pulse list (M)
CaPulseList = 1.1 * [0.1 0.5 0.75 1.0 2 3 4 5 7.5 10] * 1e-6;

%% options
opts = odeset('RelTol',1e-12,'AbsTol',1e-10);
colors = flipud(hsv(length(CaPulseList)));

%% WT figure axes

figure(1); clf; set(gcf,'Position',[20,120,1250,600],'Color','w');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
% sgtitle(['WT — ' condLabel], 'FontSize', 13, 'FontWeight', 'bold');

nexttile; axis off; % blank

axWT_P  = nexttile; hold on; xlabel('Time (s)'); ylabel('P_{Open}');
xlim([650 1000]); ylim([0 1]); set(gca,'FontSize',14); box off;

axWT_Cx = nexttile; hold on; xlabel('Time (s)'); ylabel('[Ca^{2+}]_m (\muM)');
xlim([650 1000]); ylim([0.05 0.48]); yticks(0:0.1:0.5); set(gca,'FontSize',14); box off;

axWT_Ce = nexttile; hold on; xlabel('Time (s)'); ylabel('[Ca^{2+}]_c (\muM)');
xlim([650 1000]); ylim([0 9.3]); set(gca,'FontSize',14); box off;

axWT_J  = nexttile; hold on; xlabel('Time (s)'); ylabel('J_{MCU} (nmol mg^{-1} s^{-1})');
xlim([650 1000]); set(gca,'FontSize',14); box off;

axWT_BP = nexttile; hold on; xlabel('Time (s)'); ylabel('[BPCa]_m (\muM)');
xlim([650 1000]); ylim([500 4075]); set(gca,'FontSize',14); box off;

%% Make KO figure axes

figure(2); clf; set(gcf,'Position',[500,120,1250,600],'Color','w');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
% sgtitle(['KO — ' condLabel], 'FontSize', 13, 'FontWeight', 'bold');

nexttile; axis off; % blank

axKO_P  = nexttile; hold on; xlabel('Time (s)'); ylabel('P_{Open}');
xlim([650 1000]); ylim([0 1.25]); set(gca,'FontSize',14); box off;

axKO_Cx = nexttile; hold on; xlabel('Time (s)'); ylabel('[Ca^{2+}]_m (\muM)');
xlim([650 1000]); ylim([0.05 0.6]); yticks(0:0.1:0.5); set(gca,'FontSize',14); box off;

axKO_Ce = nexttile; hold on; xlabel('Time (s)'); ylabel('[Ca^{2+}]_c (\muM)');
xlim([650 1000]); ylim([0 8]); set(gca,'FontSize',14); box off;

axKO_J  = nexttile; hold on; xlabel('Time (s)'); ylabel('J_{MCU} (nmol mg^{-1} s^{-1})');
xlim([650 1000]); ylim([0 0.065]); set(gca,'FontSize',14); box off;

axKO_BP = nexttile; hold on; xlabel('Time (s)'); ylabel('[BPCa]_m (\muM)');
xlim([650 1000]); ylim([500 4800]); set(gca,'FontSize',14); box off;

%% Ca pulses
for k = 1:length(CaPulseList)

    % Ca pulse amplitude (global)
    g_Ca_amplitude = CaPulseList(k);

    % label
    labelAmp = sprintf('%.2f \\muM', g_Ca_amplitude*1e6);

    % color
    c = colors(k,:);

    %  WT SIMULATION 
    USE_G_UNITY = 0; 
    g0  = WT_g0;
    g1  = WT_g1;
    g12 = WT_g12;

    [tWT, yWT] = run_two_parts(IC0, t_Ca_add, t_end, opts);

    Ca_e_WT = yWT(:,5);
    Ca_x_WT = yWT(:,6);
    x0_WT   = yWT(:,7);
    x1_WT   = yWT(:,8);
    x12_WT  = 1 - x0_WT - x1_WT;
    
    % Extract inactivation states
    S_O_WT = yWT(:,9);
    S_I_WT = yWT(:,10);
    S_R_WT = yWT(:,11);
    if USE_INACTIVATION
        P_open_x_WT = S_O_WT + S_R_WT;  % inactivation ON
    else
        P_open_x_WT = ones(size(S_O_WT)); % inactivation OFF
    end
    
    % Base gating
    Popen_e = g0*x0_WT + g1*x1_WT + g12*x12_WT;
    
    % Effective open probability (MICU1/2 × EMRE)
    Popen_WT = Popen_e .* P_open_x_WT;

    [~, BPCa_WT, ~] = Ca_Buffering_in(Ca_x_WT);

    Jmcu_WT = make_Jmcu(Ca_x_WT, Ca_e_WT, Popen_WT);

    plot(axWT_P,  tWT, Popen_WT, '-', 'LineWidth',1.4, 'Color',c, 'DisplayName',['WT  ' labelAmp '  (' condLabel ')']);
    plot(axWT_Cx, tWT, Ca_x_WT*1e6, '-', 'LineWidth',1.6, 'Color',c);
    plot(axWT_Ce, tWT, Ca_e_WT*1e6, '-', 'LineWidth',1.6, 'Color',c);
    plot(axWT_J,  tWT, Jmcu_WT*1e3, '-', 'LineWidth',1.4, 'Color',c);
    plot(axWT_BP, tWT, BPCa_WT*1e6, '-', 'LineWidth',1.4, 'Color',c);

    %  KO SIMULATION 
    USE_G_UNITY = 1;  
    g0  = KO_g0;
    g1  = KO_g1;
    g12 = KO_g12;

    [tKO, yKO] = run_two_parts(IC0, t_Ca_add, t_end, opts);

    Ca_e_KO = yKO(:,5);
    Ca_x_KO = yKO(:,6);
    x0_KO   = yKO(:,7);
    x1_KO   = yKO(:,8);
    x12_KO  = 1 - x0_KO - x1_KO;
    
    % For KO: constant open, so P_open is just 0.95
    Popen_KO = 1 * ones(size(x0_KO));

    [~, BPCa_KO, ~] = Ca_Buffering_in(Ca_x_KO);

    Jmcu_KO = make_Jmcu(Ca_x_KO, Ca_e_KO, Popen_KO);

    plot(axKO_P,  tKO, Popen_KO, '-', 'LineWidth',1.4, 'Color',c, 'DisplayName',['KO  ' labelAmp '  (' condLabel ')']);
    plot(axKO_Cx, tKO, Ca_x_KO*1e6, '-', 'LineWidth',1.6, 'Color',c);
    plot(axKO_Ce, tKO, Ca_e_KO*1e6, '-', 'LineWidth',1.6, 'Color',c);
    plot(axKO_J,  tKO, Jmcu_KO*1e3, '-', 'LineWidth',1.4, 'Color',c);
    plot(axKO_BP, tKO, BPCa_KO*1e6, '-', 'LineWidth',1.4, 'Color',c);

end

% =========

function [tAll, yAll] = run_two_parts(IC0, tAdd, tEnd, opts)

    global t_Ca_add
    t_Ca_add = tAdd;

    tspan1 = 0:0.1:tAdd;
    [t1, y1] = ode15s(@Cation_Dynamics_ODEs, tspan1, IC0, opts);

    IC_after = y1(end,:);
    tspan2 = tAdd:0.1:tEnd;
    [t2, y2] = ode15s(@Cation_Dynamics_ODEs, tspan2, IC_after, opts);

    tAll = [t1; t2(2:end)];
    yAll = [y1; y2(2:end,:)];
end

function J = make_Jmcu(Ca_x, Ca_e, Popen)
% J_MCU over time using your MCU function
    Mg_e   = 0.5e-3;
    dPsi   = 180e-3;

    J = (0.02).* Popen(:) .* MCU(Ca_x(:), Ca_e(:), Mg_e, dPsi);
end