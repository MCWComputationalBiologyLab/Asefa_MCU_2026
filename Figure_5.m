close all; clear; clc;
% =========================================================================
% Figure 5
%     black = MICU gating ON  (mode 'WT')
%     red   = gate forced open (mode 'FORCED_OPEN')
% =========================================================================

%% Paths:
addpath('.\Cation_Fluxes')

%% Putative biphasic [Ca2+]m regulation of MCU Gating 
USE_INACTIVATION = 0;    % 0 -> P_open_x = 1 ; 1 -> P_open_x = S_O + S_R
labelBlack = 'With Biphasic [Ca^{2+}]_m Regulation of MCU Kinetics';
labelRed   = 'No [Ca^{2+}]_m Regulation of MCU Kinetics';

%% Base parameters
p = cation_params();            % baseline

p.use_inactivation = USE_INACTIVATION;
p.t_Ca_add = 700;

options = odeset('RelTol',1e-12,'AbsTol',1e-10,'MaxStep',0.2,'NonNegative',1:6,'Refine',1);

%% Initial conditions (13-state)
IC0 = p.IC0;

t_end = 1100;
Ca_Pulse = [0.3, 0.7, 1.0, 1.5]*1e-6;

figure(1); clf; set(gcf,'Position',[20,120,1500,570],'Color','w'); hold on
axBlank = nexttile; axis(axBlank,'off');

axPopen = nexttile; hold(axPopen,'on');
ylabel(axPopen,'P_{Open} (unitless)'); xlabel(axPopen,'Time (s)');
xlim(axPopen,[650 1000]); xticks(axPopen,700:100:1100); set(axPopen,'FontSize',16);
ylim(axPopen,[0 0.3]); yticks(axPopen,0:0.05:0.15);

axCx = nexttile; hold(axCx,'on');
ylabel(axCx,'[Ca^{2+}]_m (\muM)'); xlabel(axCx,'Time (s)');
xlim(axCx,[650 1000]); xticks(axCx,700:100:1100); set(axCx,'FontSize',16);
ylim(axCx,[0.04 .12]); yticks(axCx,0.0:0.04:0.25);

axCe = nexttile; hold(axCe,'on');
ylabel(axCe,'[Ca^{2+}]_c (nM)'); xlabel(axCe,'Time (s)');
xlim(axCe,[650 1000]); xticks(axCe,700:100:1100); set(axCe,'FontSize',16);
ylim(axCe,[0 2000]); yticks(axCe,0:500:2000);

axJmcu = nexttile; hold(axJmcu,'on');
ylabel(axJmcu,'J_{MCUx} (nmol mg^{-1} s^{-1})'); xlabel(axJmcu,'Time (s)');
xlim(axJmcu,[650 1100]); set(axJmcu,'FontSize',16);

axBPCa2 = nexttile; hold(axBPCa2,'on');
ylabel(axBPCa2,'[BPCa]_m (\muM)'); xlabel(axBPCa2,'Time (s)');
xlim(axBPCa2,[650 1000]); xticks(axBPCa2,700:100:1100); set(axBPCa2,'FontSize',16);
ylim(axBPCa2,[370 1200]); yticks(axBPCa2,00:250:2000);

run_and_plot(IC0, t_end, Ca_Pulse, p, options, ...
    axPopen, axCx, axCe, axJmcu, axBPCa2, 'nM', labelBlack, labelRed);

%% =======================================================================
%  RUN 2:
%% =======================================================================

t_end_F6 = 1200;
Ca_amplitudes_F6 = [1.5 2.5 5.0 7.5 10]*1e-6;

figure(2); clf; set(gcf,'Position',[500,120,1500,570],'Color','w'); hold on
axBlank = nexttile; axis(axBlank,'off');

axPopen = nexttile; hold(axPopen,'on');
ylabel(axPopen,'P_{Open} (unitless)'); xlabel(axPopen,'Time (s)');
xlim(axPopen,[650 1000]); xticks(axPopen,700:100:1100);
ylim(axPopen,[0 1.03]); yticks(axPopen,0:0.2:1); set(axPopen,'FontSize',16);

axCx = nexttile; hold(axCx,'on');
ylabel(axCx,'[Ca^{2+}]_m (\muM)'); xlabel(axCx,'Time (s)');
xlim(axCx,[650 1000]); xticks(axCx,700:100:1100);
ylim(axCx,[0.05 0.7]); yticks(axCx,0:0.2:1.2); set(axCx,'FontSize',16);

axCe = nexttile; hold(axCe,'on');
ylabel(axCe,'[Ca^{2+}]_c (\muM)'); xlabel(axCe,'Time (s)');
xlim(axCe,[650 1000]); xticks(axCe,700:100:1100);
yticks(axCe,0:2:8); set(axCe,'FontSize',16);

axJmcu = nexttile; hold(axJmcu,'on');
ylabel(axJmcu,'J_{MCUx} (nmol mg^{-1} s^{-1})'); xlabel(axJmcu,'Time (s)');
xlim(axJmcu,[650 1000]); set(axJmcu,'FontSize',16);

axBPCa2 = nexttile; hold(axBPCa2,'on');
ylabel(axBPCa2,'[BPCa]_m (\muM)'); xlabel(axBPCa2,'Time (s)');
xlim(axBPCa2,[650 1000]); xticks(axBPCa2,700:100:1100); set(axBPCa2,'FontSize',16);
ylim(axBPCa2,[500 5400]); yticks(axBPCa2,500:1500:5000); set(axBPCa2,'FontSize',16);

run_and_plot(IC0, t_end_F6, Ca_amplitudes_F6, p, options, ...
    axPopen, axCx, axCe, axJmcu, axBPCa2, 'uM', labelBlack, labelRed);

%%  local simul
function run_and_plot(IC, t_end, Ca_amps, p0, options, ...
                      axPopen, axCx, axCe, axJmcu, axBPCa2, ceUnit, labelBlack, labelRed)
    for i = 1:numel(Ca_amps)
        amp = Ca_amps(i);

        % (A) MICU gating ON (black) 
        pM = p0; pM.mode = 'WT'; pM.Ca_Pulse = amp;
        odeM = @(t,y) cation_dynamics_odes(t,y,pM);
        [tM, yM] = ode15s(odeM, [0 t_end], IC, options);

        Ca_e_M = yM(:,5); Ca_x_M = yM(:,6);
        x0=yM(:,7); x1=yM(:,8); x12=1-x0-x1;
        S_O=yM(:,11); S_R=yM(:,13);
        if pM.use_inactivation, PoxM=S_O+S_R; else, PoxM=ones(size(S_O)); end
        Popen_M = (pM.g0*x0 + pM.g1*x1 + pM.g12*x12) .* PoxM;   % effective gate
        fM = cation_fluxes(tM, yM, pM);
        [~, BPCaM] = Ca_Buffering_in(Ca_x_M);

        vis = 'off'; if i==1, vis = 'on'; end
        plot(axPopen, tM, Popen_M, '-','Color','k','LineWidth',1.1, ...
             'DisplayName',labelBlack,'HandleVisibility',vis);
        plot(axCx,   tM, Ca_x_M*1e6,      '-','Color','k','LineWidth',1.3);
        plot(axCe,   tM, ce_scale(Ca_e_M,ceUnit), '-','Color','k','LineWidth',1.1);
        plot(axJmcu, tM, fM.J_MCU*1e3,    '-','Color','k','LineWidth',1.0);
        plot(axBPCa2,tM, BPCaM*1e6,       '-','Color','k','LineWidth',1.0);

        % (B) Gate forced open (red)
        pU = p0; pU.mode = 'FORCED_OPEN'; pU.Ca_Pulse = amp;
        odeU = @(t,y) cation_dynamics_odes(t,y,pU);
        [tU, yU] = ode15s(odeU, [0 t_end], IC, options);

        Ca_e_U = yU(:,5); Ca_x_U = yU(:,6);
        S_O=yU(:,11); S_R=yU(:,13);
        if pU.use_inactivation, PoxU=S_O+S_R; else, PoxU=ones(size(S_O)); end
        Popen_U = PoxU;    % forced-open panel
        fU = cation_fluxes(tU, yU, pU);
        [~, BPCaU] = Ca_Buffering_in(Ca_x_U);

        plot(axPopen, tU, Popen_U, '-','Color','r','LineWidth',1.1, ...
             'DisplayName',labelRed,'HandleVisibility',vis);
        plot(axCx,   tU, Ca_x_U*1e6,      '-','Color','r','LineWidth',1.1);
        plot(axCe,   tU, ce_scale(Ca_e_U,ceUnit), '-','Color','r','LineWidth',1.1);
        plot(axJmcu, tU, fU.J_MCU*1e3,    '-','Color','r','LineWidth',1.0);
        plot(axBPCa2,tU, BPCaU*1e6,       '-','Color','r','LineWidth',1.0);
    end
end

function v = ce_scale(Ca_e, unit)
    if strcmpi(unit,'nM'), v = Ca_e*1e9; else, v = Ca_e*1e6; end
end
