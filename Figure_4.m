close all; clear; clc;
% =========================================================================
% Figure 4 
%   WT vs KO(forced-open) MCU response to a range of Ca2+ pulses.
% =========================================================================

%% Paths:
addpath('.\Cation_Fluxes')

%%  Putative biphasic [Ca2+]m regulation of MCU Gating
% 0 -> P_open_x = 1 (inactivation OFF) ; 1 -> P_open_x = S_O + S_R (ON)
USE_INACTIVATION = 1;
if USE_INACTIVATION
    condLabel = 'With Biphasic [Ca^{2+}]_m Regulation of MCU Kinetics';
else
    condLabel = 'No [Ca^{2+}]_m Regulation of MCU Kinetics';
end

%% Base parameters
p = cation_params();          % baseline == original Cation_Dynamics_ODEs.m

p.use_inactivation = USE_INACTIVATION;
t_end    = 1100;
t_Ca_add = 700;

%% Initial conditions (13-state)
IC0 = p.IC0;

%% Ca pulse list (M)
CaPulseList = [0.1 0.5 0.75 1.0 2 3 4 5 7.5 9.5] * 1e-6;

opts   = odeset('RelTol',1e-12,'AbsTol',1e-10);
colors = flipud(hsv(numel(CaPulseList)));

%% WT figure axes
figure(1); clf; set(gcf,'Position',[20,120,1250,600],'Color','w');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

nexttile; axis off;
axWT_P  = nexttile; hold on; xlabel('Time (s)'); ylabel('P_{Open}');
xlim([650 1000]); ylim([0 1]); set(gca,'FontSize',14); box off;
axWT_Cx = nexttile; hold on; xlabel('Time (s)'); ylabel('[Ca^{2+}]_m (\muM)');
xlim([650 1000]); ylim([0.05 0.6]); yticks(0:0.2:0.6); set(gca,'FontSize',14); box off;
axWT_Ce = nexttile; hold on; xlabel('Time (s)'); ylabel('[Ca^{2+}]_c (\muM)');
xlim([650 1000]); ylim([0 9.3]); set(gca,'FontSize',14); box off;
axWT_J  = nexttile; hold on; xlabel('Time (s)'); ylabel('J_{MCUx} (nmol mg^{-1} s^{-1})');
xlim([650 1000]);  yticks(0:0.02:0.06); set(gca,'FontSize',14); box off;
axWT_BP = nexttile; hold on; xlabel('Time (s)'); ylabel('[BPCa]_m (\muM)');
xlim([650 1000]); ylim([500 4700]); set(gca,'FontSize',14); box off;

%% KO figure axes 
figure(2); clf; set(gcf,'Position',[500,120,1250,600],'Color','w');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

nexttile; axis off;
axKO_P  = nexttile; hold on; xlabel('Time (s)'); ylabel('P_{Open}');
xlim([650 1000]); ylim([0 1.2]); yticks(0:0.5:1);set(gca,'FontSize',14); box off;
axKO_Cx = nexttile; hold on; xlabel('Time (s)'); ylabel('[Ca^{2+}]_m (\muM)');
xlim([650 1000]); ylim([0.05 0.65]); yticks(0:0.2:0.6); set(gca,'FontSize',14); box off;
axKO_Ce = nexttile; hold on; xlabel('Time (s)'); ylabel('[Ca^{2+}]_c (\muM)');
xlim([650 1000]); ylim([0 8]); set(gca,'FontSize',14); box off;
axKO_J  = nexttile; hold on; xlabel('Time (s)'); ylabel('J_{MCUx} (nmol mg^{-1} s^{-1})');
xlim([650 1000]); ylim([0 0.065]); yticks(0:0.02:0.06); set(gca,'FontSize',14); box off;
axKO_BP = nexttile; hold on; xlabel('Time (s)'); ylabel('[BPCa]_m (\muM)');
xlim([650 1000]); ylim([500 5200]); set(gca,'FontSize',14); box off;

%%  Sweep Ca pulses 
for k = 1:numel(CaPulseList)
    amp = CaPulseList(k);
    labelAmp = sprintf('%.2f \\muM', amp*1e6);
    c = colors(k,:);

    %  WT simulation 
    pWT = p; pWT.mode = 'WT'; pWT.Ca_Pulse = amp;
    [tWT, yWT] = run_two_parts(IC0, t_Ca_add, t_end, pWT, opts);

    Ca_e_WT = yWT(:,5);  Ca_x_WT = yWT(:,6);
    x0 = yWT(:,7); x1 = yWT(:,8); x12 = 1 - x0 - x1;
    S_O = yWT(:,11); S_R = yWT(:,13);
    if USE_INACTIVATION, Pox = S_O + S_R; else, Pox = ones(size(S_O)); end
    Popen_WT = (pWT.g0*x0 + pWT.g1*x1 + pWT.g12*x12) .* Pox;   % effective gate
    fWT = cation_fluxes(tWT, yWT, pWT);
    [~, BPCa_WT] = Ca_Buffering_in(Ca_x_WT);

    plot(axWT_P,  tWT, Popen_WT,       '-','LineWidth',1.4,'Color',c, ...
         'DisplayName',['WT  ' labelAmp '  (' condLabel ')']);
    plot(axWT_Cx, tWT, Ca_x_WT*1e6,    '-','LineWidth',1.6,'Color',c);
    plot(axWT_Ce, tWT, Ca_e_WT*1e6,    '-','LineWidth',1.6,'Color',c);
    plot(axWT_J,  tWT, fWT.J_MCU*1e3,  '-','LineWidth',1.4,'Color',c);
    plot(axWT_BP, tWT, BPCa_WT*1e6,    '-','LineWidth',1.4,'Color',c);

    % ----- KO simulation (forced-open gate) -----
    pKO = p; pKO.mode = 'FORCED_OPEN'; pKO.Ca_Pulse = amp;
    [tKO, yKO] = run_two_parts(IC0, t_Ca_add, t_end, pKO, opts);

    Ca_e_KO = yKO(:,5);  Ca_x_KO = yKO(:,6);
    S_O = yKO(:,11); S_R = yKO(:,13);
    if USE_INACTIVATION, PoxK = S_O + S_R; else, PoxK = ones(size(S_O)); end
    Popen_KO = PoxK;    % KO weights are unity -> MICU gate = 1 (panel convention)
    fKO = cation_fluxes(tKO, yKO, pKO);
    [~, BPCa_KO] = Ca_Buffering_in(Ca_x_KO);

    plot(axKO_P,  tKO, Popen_KO,       '-','LineWidth',1.4,'Color',c, ...
         'DisplayName',['KO  ' labelAmp '  (' condLabel ')']);
    plot(axKO_Cx, tKO, Ca_x_KO*1e6,    '-','LineWidth',1.6,'Color',c);
    plot(axKO_Ce, tKO, Ca_e_KO*1e6,    '-','LineWidth',1.6,'Color',c);
    plot(axKO_J,  tKO, fKO.J_MCU*1e3,  '-','LineWidth',1.4,'Color',c);
    plot(axKO_BP, tKO, BPCa_KO*1e6,    '-','LineWidth',1.4,'Color',c);
end

%% local simulation  
function [tAll, yAll] = run_two_parts(IC0, tAdd, tEnd, p, opts)
    p.t_Ca_add = tAdd;
    ode = @(t,y) cation_dynamics_odes(t, y, p);
    [t1, y1] = ode15s(ode, 0:0.1:tAdd,    IC0,        opts);
    [t2, y2] = ode15s(ode, tAdd:0.1:tEnd, y1(end,:),  opts);
    tAll = [t1; t2(2:end)];
    yAll = [y1; y2(2:end,:)];
end
