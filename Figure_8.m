close all; clear; clc;
% =========================================================================
% Figure 8  (refactored)  -- 3x3 dashboard, produced TWICE:
%    Figure 1 -> USE_INACTIVATION = 0   (P_open = 1)
%    Figure 2 -> USE_INACTIVATION = 1   (P_open = S_O + S_R)
% Cytosolic [Mg2+] conditions; time-course + Ca dose-response.

% =========================================================================

%% Paths:
addpath('.\Cation_Fluxes')

%%  Base parameters (WT)
p = cation_params();

%%  Time-course parameters
t_Ca_add_TC = 500;   t_end_TC = 800;   Ca_amp_TC = 20e-6;

%%  Dose-response parameters 
CaPulseList = [0 0.25 0.5 1 2 3 4 5 6 8 10 15 20]*1e-6;
n_pulses    = numel(CaPulseList);
Ca_pulse_uM = CaPulseList*1e6;
t_Ca_add_DR = 50;    t_end_DR = 300;   t_ss = 290;

%% ---- 4 cytosolic Mg conditions
Mg_levels_mM = [0, 1, 5, 10];
n_Mg = numel(Mg_levels_mM);
Mg_colors = { [0.50 0.00 0.60], [1.00 0.60 0.70], [0.90 0.10 0.50], [1.00 0.55 0.00] };

%% Initial conditions (13-state)
IC0 = p.IC0;

opts = odeset('RelTol',1e-12,'AbsTol',1e-10,'NonNegative',1:13);

%%  Inactivation conditions
inact_vals   = [0, 1];
inact_titles = { 'Inactivation OFF  -  P_{open} = 1', ...
                 'Inactivation ON   -  P_{open} = S_O + S_R' };
n_inact = numel(inact_vals);

%% ========================================================================
for iInact = 1:n_inact
    p.use_inactivation = inact_vals(iInact);

    Ca_e_ss = zeros(n_pulses, n_Mg);
    Ca_x_ss = zeros(n_pulses, n_Mg);

    fh = figure(iInact); clf(fh);
    set(fh,'Position',[30+(iInact-1)*100, 60, 1350, 900],'Color','w');
    sgtitle(fh, inact_titles{iInact}, 'FontSize',17,'FontWeight','bold');
    tl = tiledlayout(fh,3,3,'Padding','compact','TileSpacing','compact');

    ax1 = nexttile(tl,1);  axis(ax1,'off');
    ax2 = nexttile(tl,2);  hold(ax2,'on'); box(ax2,'off');
    ylabel(ax2,'[Ca^{2+}]_c (\muM)'); xlabel(ax2,'Time (s)');
    xlim(ax2,[450,800]); ylim(ax2,[-0.5,21]); yticks(ax2,0:5:20); set(ax2,'FontSize',16);
    ax3 = nexttile(tl,3);  hold(ax3,'on'); box(ax3,'off');
    ylabel(ax3,'P_{open}'); xlabel(ax3,'Time (s)');
    xlim(ax3,[450,800]); ylim(ax3,[0,1.05]); yticks(ax3,0:0.25:1); set(ax3,'FontSize',16);
    ax4 = nexttile(tl,4);  hold(ax4,'on'); box(ax4,'off');
    ylabel(ax4,'J_{MCUx} (nmol/mg/s)'); xlabel(ax4,'Time (s)');
    xlim(ax4,[450,800]); ylim(ax4,[0,0.105]); yticks(ax4,0:0.025:0.1); set(ax4,'FontSize',16);
    ax5 = nexttile(tl,5);  hold(ax5,'on'); box(ax5,'off');
    ylabel(ax5,'[Ca^{2+}]_m (\muM)'); xlabel(ax5,'Time (s)');
    xlim(ax5,[450,800]); ylim(ax5,[0,2]); yticks(ax5,0:0.5:2); set(ax5,'FontSize',16);
    ax6 = nexttile(tl,6);  hold(ax6,'on'); box(ax6,'off');
    ylabel(ax6,'[CaBP] (mM)'); xlabel(ax6,'Time (s)');
    xlim(ax6,[450,800]); ylim(ax6,[0,11]); yticks(ax6,0:2.5:10); set(ax6,'FontSize',16);
    ax7 = nexttile(tl,7);  hold(ax7,'on'); box(ax7,'off');
    xlabel(ax7,'Ca^{2+} Pulse (\muM)'); ylabel(ax7,'SS [Ca^{2+}]_c (\muM)'); set(ax7,'FontSize',16);
    ax8 = nexttile(tl,8);  hold(ax8,'on'); box(ax8,'off');
    xlabel(ax8,'Ca^{2+} Pulse (\muM)'); ylabel(ax8,'SS [Ca^{2+}]_m (\muM)'); yticks(ax8,0:0.5:2); set(ax8,'FontSize',16);
    ax9 = nexttile(tl,9);  hold(ax9,'on'); box(ax9,'off');
    xlabel(ax9,'SS [Ca^{2+}]_c (\muM)'); ylabel(ax9,'SS [Ca^{2+}]_m (\muM)');  xticks(ax9,0:2:10); yticks(ax9,0:0.5:2); set(ax9,'FontSize',16);

    for m = 1:n_Mg
        p = p;
        p.Mg_e = Mg_levels_mM(m)*1e-3;
        Mg_nom = Mg_levels_mM(m);
        lbl = sprintf('[Mg^{2+}]_c = %.4g mM', Mg_nom);

        %% (A) Time-course
        p.Ca_Pulse = Ca_amp_TC;
        [tTC, yTC] = run_two_parts(IC0, t_Ca_add_TC, t_end_TC, p, opts);

        Ca_e_tc = yTC(:,5);  Ca_x_tc = yTC(:,6);
        fTC = cation_fluxes(tTC, yTC, p);        % consistent J_MCU + effective gate
        P_open_tc = fTC.G;
        [~, CaBP_tc] = Ca_Buffering_in(Ca_x_tc);

        plot(ax2, tTC, Ca_e_tc*1e6,  '-','LineWidth',2.1,'Color',Mg_colors{m},'DisplayName',lbl);
        plot(ax3, tTC, P_open_tc,    '-','LineWidth',2.1,'Color',Mg_colors{m},'HandleVisibility','off');
        plot(ax4, tTC, fTC.J_MCU*1e3,'-','LineWidth',2.1,'Color',Mg_colors{m},'HandleVisibility','off');
        plot(ax5, tTC, Ca_x_tc*1e6,  '-','LineWidth',2.1,'Color',Mg_colors{m},'HandleVisibility','off');
        plot(ax6, tTC, CaBP_tc*1e3,  '-','LineWidth',2.1,'Color',Mg_colors{m},'HandleVisibility','off');

        %% (B) Ca dose-response sweep
        for k = 1:n_pulses
            p.Ca_Pulse = CaPulseList(k);
            [tDR, yDR]   = run_two_parts(IC0, t_Ca_add_DR, t_end_DR, p, opts);
            [~, idx_ss]  = min(abs(tDR - t_ss));
            Ca_e_ss(k,m) = yDR(idx_ss,5);
            Ca_x_ss(k,m) = yDR(idx_ss,6);
        end

        plot(ax7, Ca_pulse_uM, Ca_e_ss(:,m)*1e6, 'o-','LineWidth',2.1,'MarkerSize',2.5, ...
             'Color',Mg_colors{m},'MarkerFaceColor',Mg_colors{m},'HandleVisibility','off');
        plot(ax8, Ca_pulse_uM, Ca_x_ss(:,m)*1e6, 'o-','LineWidth',2.1,'MarkerSize',2.5, ...
             'Color',Mg_colors{m},'MarkerFaceColor',Mg_colors{m},'HandleVisibility','off');
        plot(ax9, Ca_e_ss(:,m)*1e6, Ca_x_ss(:,m)*1e6, 'o-','LineWidth',2.1,'MarkerSize',2.5, ...
             'Color',Mg_colors{m},'MarkerFaceColor',Mg_colors{m},'HandleVisibility','off');
    end

    hLines = flipud(findobj(ax2,'Type','Line'));
    legend(ax1, hLines, {hLines.DisplayName}, ...
           'Orientation','vertical','Location','best','FontSize',14,'Box','off');
    ax1.Visible = 'off';
    linkaxes([ax2 ax3 ax4 ax5 ax6],'x');
end

%% local simulation
function [tAll, yAll] = run_two_parts(IC0, tAdd, tEnd, p, opts)
    p.t_Ca_add = tAdd;
    ode = @(t,y) cation_dynamics_odes(t,y,p);
    [t1, y1] = ode15s(ode, 0:0.1:tAdd,    IC0,       opts);
    [t2, y2] = ode15s(ode, tAdd:0.1:tEnd, y1(end,:), opts);
    tAll = [t1; t2(2:end)];
    yAll = [y1; y2(2:end,:)];
end
