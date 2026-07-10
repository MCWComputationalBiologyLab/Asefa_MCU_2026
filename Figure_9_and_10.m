close all; clear; clc;
% =========================================================================
% Figure 9 & 10
%   WT Ca2+ dose-response, two conditions:
%     cond 1 : No  [Ca2+]m regulation      (p.use_inactivation = 0, P_open_x = 1)
%     cond 2 : With biphasic [Ca2+]m reg.  (p.use_inactivation = 1, P_open_x = S_O+S_R)
% =========================================================================

%% Paths:
addpath('.\Cation_Fluxes')

%% Base parameters
p = cation_params();          % baseline

% --- Pulse settings ---
p.t_Ca_add = 50;

%% Simulation parameters
t_end    = 300;
t_Ca_add = 50;
t_ss     = 290;

%% Initial conditions (13-state)
IC0 = p.IC0;

%% Ca pulse list (M)
CaPulseList = [0 .1 .25 .5 .75 1 2 3 4 5 6 7 8 9 10 12.5 15 17.5 20] * 1e-6;
% CaPulseList = [0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1 1.25 1.5 1.75 2 2.5 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20] * 1e-6;
n_pulses    = length(CaPulseList);
Ca_pulse_uM = CaPulseList * 1e6;

%  FIGURE 1 (time-courses) 
%   Col 1 = [Ca2+]_c    Col 2 = P_Open    Col 3 = J_MCUx    Col 4 = [Ca2+]_m

Fig1_XLim   = [0 200];                      % [] -> auto
Fig1_XTicks = [0 50 100 150 200];           % [] -> auto

%  Y axis — one entry per column:
Fig1_YLim{1}   = [0 20];                 Fig1_YTicks{1} = [0 5 10 15 20];
Fig1_YLim{2}   = [0 1];                  Fig1_YTicks{2} = [0 0.25 0.5 0.75 1.0];
Fig1_YLim{3}   = [0 0.1];                Fig1_YTicks{3} = [];
Fig1_YLim{4}   = [0 1.61];               Fig1_YTicks{4} = (0:0.5:2);

%  FIGURE 2 
%  12 panels, row-major tile order (1..12):
%   1=A(Peak Ca_c)   2=B(Peak P_open)   3=C(Peak J_MCU)   4=D(SS Ca_m)
%   5=E(SS Ca_c)     6=F(SS P_open)     7=G(SS J_MCU)     8=H(SS Ca_m vs SS Ca_c)
%   9=I(AUC Ca_c)   10=J(AUC P_open)  11=K(AUC J_MCU)   12=L(AUC Ca_m vs pulse)
%              tile:   1         2          3           4
Fig2_XLim  = {   [],       [],        [],         [],   ...
                 [],       [],        [],         [],   ...
                 [],       [],        [],         []   };
Fig2_XTick = { (0:5:20), (0:5:20),  (0:5:20),   (0:5:20), ...
               (0:5:20), (0:5:20),  (0:5:20),   (0:.5:2), ...
               (0:5:20), (0:5:20),  (0:5:20),   (0:5:20) };


Fig2_YLim  = {   [],       [0 1],     [],         [0 1.63], ...
                 [],       [0 0.165],  [],         [0 1.62], ...
                 [0 1370],  [],        [0 4.6],      [0 340] };

Fig2_YTick = { (0:5:20),    (0:0.25:1), (0:0.025:.1), (0:0.5:2), ...
               (0:0.5:2),   (0:0.05:.15), [],         (0:0.5:2), ...
               (0:300:1200),(0:25:100), (0:4),        (0:100:400) };

%% Solver & colours
opts   = odeset('RelTol', 1e-12, 'AbsTol', 1e-10);
colors = flipud(hsv(n_pulses));

%% Bright-to-faded ramps for panel I (phase-plane)
fade = linspace(0.75, 0, n_pulses)';
colors_PP = zeros(2, n_pulses, 3);
for k = 1:n_pulses
    colors_PP(1, k, :) = [fade(k), fade(k), 1       ];   % blue ramp
    colors_PP(2, k, :) = [1,        fade(k), fade(k)];   % red ramp
end

%% Condition metadata
condLabels   = {'Without Matrix Ca^{2+} Regulation', 'With Matrix Ca^{2+} Regulation'};
condColors   = {'b', 'r'};
lineStyles   = {'-', '-'};
markerStyles = {'o', 's'};

%% ========================================================================
%  FIGURE 1:
%% ========================================================================
figure(1); clf;
set(gcf, 'Position', [30, 80, 1500, 700], 'Color', 'w');

tl1_outer = tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

rowTitles = {'Without Matrix Ca^{2+} Regulation of MCUx Gating', 'With Matrix Ca^{2+} Regulation of MCUx Gating'};
ylabels1  = {'[Ca^{2+}]_c (\muM)', 'P_{Open}', ...
             'J_{MCUx} (nmo/mg/s)', '[Ca^{2+}]_m (\muM)'};

ax1 = gobjects(2, 4);
tl1_inner = gobjects(1, 2);

for row = 1:2
    tl1_inner(row) = tiledlayout(tl1_outer, 1, 4, ...
        'Padding', 'compact', 'TileSpacing', 'compact');
    tl1_inner(row).Layout.Tile = row;
    title(tl1_inner(row), rowTitles{row}, 'FontSize', 18, 'FontWeight', 'bold');

    for col = 1:4
        ax1(row, col) = nexttile(tl1_inner(row));
        hold on; box on;
        xlabel('Time (s)', 'FontSize', 14);
        ylabel(ylabels1{col}, 'FontSize', 14);
        if ~isempty(Fig1_XLim),      xlim(ax1(row,col), Fig1_XLim);      end
        if ~isempty(Fig1_YLim{col}), ylim(ax1(row,col), Fig1_YLim{col}); end
        set(gca, 'FontSize', 16);
    end
end

%% ========================================================================
%  FIGURE 2: 
%% ========================================================================
figure(2); clf;
set(gcf, 'Position', [30, 30, 1500, 900], 'Color', 'w');
tl2 = tiledlayout(3, 4, 'Padding', 'compact', 'TileSpacing', 'compact');

ylabels2 = {'Peak [Ca^{2+}]_c (\muM)', ...
            'Peak P_{Open}', ...
            'Peak J_{MCUx} (nmo/mg/s)', ...
            'SS [Ca^{2+}]_m (\muM)', ...
            'SS [Ca^{2+}]_c (\muM)', ...
            'SS P_{Open}', ...
            'SS J_{MCUx} (nmo/mg/s)', ...
            'SS [Ca^{2+}]_m (\muM)', ...
            'AUC [Ca^{2+}]_c (\muM{\cdot}s)', ...
            'AUC P_{Open} (s)', ...
            'AUC J_{MCUx} (nmo/mg)', ...
            'AUC [Ca^{2+}]_m (\muM{\cdot}s)'};

xlabels2 = {'Ca^{2+} Pulse (\muM)', ...
            'Ca^{2+} Pulse (\muM)', ...
            'Ca^{2+} Pulse (\muM)', ...
            'Ca^{2+} Pulse (\muM)', ...
            'Ca^{2+} Pulse (\muM)', ...
            'Ca^{2+} Pulse (\muM)', ...
            'Ca^{2+} Pulse (\muM)', ...
            'SS [Ca^{2+}]_c (\muM)', ...
            'Ca^{2+} Pulse (\muM)', ...
            'Ca^{2+} Pulse (\muM)', ...
            'Ca^{2+} Pulse (\muM)', ...
            'Ca^{2+} Pulse (\muM)'};

ax2 = gobjects(1, 12);
for pIdx = 1:12
    ax2(pIdx) = nexttile(tl2);
    hold on; box on;
    xlabel(xlabels2{pIdx}, 'FontSize', 16);
    ylabel(ylabels2{pIdx}, 'FontSize', 16);
    set(gca, 'FontSize', 14);
    if ~isempty(Fig2_XLim{pIdx}),  xlim(ax2(pIdx),  Fig2_XLim{pIdx});  end
    if ~isempty(Fig2_YLim{pIdx}),  ylim(ax2(pIdx),  Fig2_YLim{pIdx});  end
    if ~isempty(Fig2_XTick{pIdx}), xticks(ax2(pIdx), Fig2_XTick{pIdx}); end
    if ~isempty(Fig2_YTick{pIdx}), yticks(ax2(pIdx), Fig2_YTick{pIdx}); end
end

%% ========================================================================
%  PREALLOCATE
%% ========================================================================
Ca_e_ss     = zeros(n_pulses, 2);
Ca_x_ss     = zeros(n_pulses, 2);
P_open_peak = zeros(n_pulses, 2);
Jmcu_peak   = zeros(n_pulses, 2);
AUC_Popen   = zeros(n_pulses, 2);
AUC_Jmcu    = zeros(n_pulses, 2);
AUC_Ca_e    = zeros(n_pulses, 2);
AUC_Ca_x    = zeros(n_pulses, 2);
Ca_e_peak   = zeros(n_pulses, 2);
Ca_x_peak   = zeros(n_pulses, 2);
P_open_ss   = zeros(n_pulses, 2);
Jmcu_ss     = zeros(n_pulses, 2);

%% ========================================================================
%  MAIN LOOP
%% ========================================================================
for cond = 1:2

    fprintf('\n=== Condition %d: %s ===\n', cond, condLabels{cond});

    % cond 1 -> no regulation (P_open_x = 1) ; cond 2 -> biphasic (S_O + S_R)
    p.use_inactivation = cond - 1;

    for k = 1:n_pulses

        amp = CaPulseList(k);
        fprintf('  Pulse %d/%d: %.2f uM\n', k, n_pulses, amp*1e6);
        c    = colors(k, :);
        c_PP = squeeze(colors_PP(cond, k, :))'; 

        %% --- RUN ODE (WT gating, shared dynamics) ---
        pWT = p; pWT.mode = 'WT'; pWT.Ca_Pulse = amp;
        [tWT, yWT] = run_two_parts(IC0, t_Ca_add, t_end, pWT, opts);

        %% --- STATES + EFFECTIVE GATE / FLUX (shared recompute) ---
        Ca_e_WT = yWT(:, 5);
        Ca_x_WT = yWT(:, 6);
        fWT      = cation_fluxes(tWT, yWT, pWT);
        Popen_WT = fWT.G;          % effective gate = (g0 x0 + g1 x1 + g12 x12) .* P_open_x
        Jmcu_WT  = fWT.J_MCU;      % to match the driving J_MCU

        %% FIGURE 1
        plot(ax1(cond,1), tWT, Ca_e_WT*1e6,  '-', 'LineWidth', 1.5, 'Color', c);
        plot(ax1(cond,2), tWT, Popen_WT,     '-', 'LineWidth', 1.4, 'Color', c);
        plot(ax1(cond,3), tWT, Jmcu_WT*1e3,  '-', 'LineWidth', 1.4, 'Color', c);
        plot(ax1(cond,4), tWT, Ca_x_WT*1e6,  '-', 'LineWidth', 1.5, 'Color', c);

        %% FIGURE 2 panel L (tile 12): AUC Ca_m vs Ca Pulse 
        % (plotted after loop in post-processing section)

        %%  STEADY-STATE & PEAKS 
        [~, idx_ss]          = min(abs(tWT - t_ss));
        Ca_e_ss(k, cond)     = Ca_e_WT(idx_ss);
        Ca_x_ss(k, cond)     = Ca_x_WT(idx_ss);
        idx_after            = tWT > t_Ca_add;
        P_open_peak(k, cond) = max(Popen_WT(idx_after));
        Ca_e_peak(k, cond)   = max(Ca_e_WT(idx_after)) * 1e6;
        Ca_x_peak(k, cond)   = max(Ca_x_WT(idx_after)) * 1e6;
        Jmcu_peak(k, cond)   = max(Jmcu_WT(idx_after)) * 1e3;

        %%  AUC 
        t_post             = tWT(idx_after);
        AUC_Popen(k, cond) = trapz(t_post, Popen_WT(idx_after));
        AUC_Jmcu(k, cond)  = trapz(t_post, Jmcu_WT(idx_after) * 1e3);
        AUC_Ca_e(k, cond)  = trapz(t_post, Ca_e_WT(idx_after) * 1e6);
        AUC_Ca_x(k, cond)  = trapz(t_post, Ca_x_WT(idx_after) * 1e6);

        %%  SS P_open & SS J_MCU 
        P_open_ss(k, cond) = Popen_WT(idx_ss);
        Jmcu_ss(k, cond)   = Jmcu_WT(idx_ss) * 1e3;

    end  % pulse loop

    cc = condColors{cond};
    ms = markerStyles{cond};

    %%  FIGURE 2:
    % Column 1: [Ca2+]_c
    plot(ax2(1), Ca_pulse_uM, Ca_e_peak(:,cond), [ms '-'], ...
         'LineWidth', 2, 'MarkerSize', 2, 'Color', cc, 'MarkerFaceColor', cc, ...
         'DisplayName', condLabels{cond});
    plot(ax2(5), Ca_pulse_uM, Ca_e_ss(:,cond)*1e6, [ms '-'], ...
         'LineWidth', 2, 'MarkerSize', 2, 'Color', cc, 'MarkerFaceColor', cc, ...
         'HandleVisibility', 'off');
    plot(ax2(9), Ca_pulse_uM, AUC_Ca_e(:,cond), [ms '-'], ...
         'LineWidth', 2, 'MarkerSize', 2, 'Color', cc, 'MarkerFaceColor', cc, ...
         'HandleVisibility', 'off');

    % Column 2: P_Open
    plot(ax2(2), Ca_pulse_uM, P_open_peak(:,cond), [ms '-'], ...
         'LineWidth', 2, 'MarkerSize', 2, 'Color', cc, 'MarkerFaceColor', cc, ...
         'HandleVisibility', 'off');
    plot(ax2(6), Ca_pulse_uM, P_open_ss(:,cond), [ms '-'], ...
         'LineWidth', 2, 'MarkerSize', 2, 'Color', cc, 'MarkerFaceColor', cc, ...
         'DisplayName', condLabels{cond});
    plot(ax2(10), Ca_pulse_uM, AUC_Popen(:,cond), [ms '-'], ...
         'LineWidth', 2, 'MarkerSize', 2, 'Color', cc, 'MarkerFaceColor', cc, ...
         'HandleVisibility', 'off');

    % Column 3: J_MCU
    plot(ax2(3), Ca_pulse_uM, Jmcu_peak(:,cond), [ms '-'], ...
         'LineWidth', 2, 'MarkerSize', 2, 'Color', cc, 'MarkerFaceColor', cc, ...
         'HandleVisibility', 'off');
    plot(ax2(7), Ca_pulse_uM, Jmcu_ss(:,cond), [ms '-'], ...
         'LineWidth', 2, 'MarkerSize', 2, 'Color', cc, 'MarkerFaceColor', cc, ...
         'HandleVisibility', 'off');
    plot(ax2(11), Ca_pulse_uM, AUC_Jmcu(:,cond), [ms '-'], ...
         'LineWidth', 2, 'MarkerSize', 2, 'Color', cc, 'MarkerFaceColor', cc, ...
         'HandleVisibility', 'off');

    % Column 4: [Ca2+]_m / extras
    plot(ax2(4), Ca_pulse_uM, Ca_x_ss(:,cond)*1e6, [ms '-'], ...
         'LineWidth', 2, 'MarkerSize', 2, 'Color', cc, 'MarkerFaceColor', cc, ...
         'HandleVisibility', 'off');
    plot(ax2(8), Ca_e_ss(:,cond)*1e6, Ca_x_ss(:,cond)*1e6, [ms '-'], ...
         'LineWidth', 2, 'MarkerSize', 2, 'Color', cc, 'MarkerFaceColor', cc, ...
         'HandleVisibility', 'off');

end  % condition loop

%%  Legend for panel A 
% legend(ax2(6), 'Location', 'southeast', 'FontSize', 15, 'Box', 'off');

%% ---- Plot panel L (tile 12): Ca Pulse vs AUC Ca_m ----
for cond = 1:2
    cc = condColors{cond};
    ms = markerStyles{cond};
    hVis = 'off';
    if cond == 1, hVis = 'on'; end
    plot(ax2(12), Ca_pulse_uM, AUC_Ca_x(:,cond), [ms '-'], ...
         'LineWidth', 2, 'MarkerSize', 2, 'Color', cc, 'MarkerFaceColor', cc, ...
         'DisplayName', condLabels{cond}, 'HandleVisibility', hVis);
end

%% ========================================================================
%  FIGURE 2 PANEL LABELS  (A-L)
%% ========================================================================
fig2_panelLabels = {'A','D','G','J','B','E','H','K','C','F','I','L'};
for pIdx = 1:12
    text(ax2(pIdx), 0.018, 1.0, fig2_panelLabels{pIdx}, ...
         'Units', 'normalized', ...
         'FontSize', 20, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'left', ...
         'VerticalAlignment', 'top');
end

%%  Optional colorbar for Figure 1 
figure(1);
panelLabels = {'A','B','C','D'; ...   % row 1: tiles (1,1)..(1,4)
               'E','F','G','H'};      % row 2: tiles (2,1)..(2,4)

for row = 1:2
    for col = 1:4
        ax = ax1(row, col);

        %  Limits 
        if ~isempty(Fig1_XLim),      xlim(ax, Fig1_XLim);      end
        if ~isempty(Fig1_YLim{col}), ylim(ax, Fig1_YLim{col}); end

        %  X ticks (shared time axis) 
        if ~isempty(Fig1_XTicks)
            set(ax, 'XTick', Fig1_XTicks);
        end

        %  Y ticks (per column) 
        if ~isempty(Fig1_YTicks{col})
            set(ax, 'YTick', Fig1_YTicks{col});
        end

        %  Panel label (upper-left, inside axes box)
        text(ax, 0.03, 0.99, panelLabels{row, col}, ...
             'Units', 'normalized', ...
             'FontSize', 20, 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'left', ...
             'VerticalAlignment', 'top');
    end
end

%% ========================================================================
%  local simulation
%% ========================================================================
function [tAll, yAll] = run_two_parts(IC0, tAdd, tEnd, p, opts)
    p.t_Ca_add = tAdd;
    ode = @(t,y) cation_dynamics_odes(t, y, p);
    [t1, y1] = ode15s(ode, 0:0.1:tAdd,    IC0,        opts);
    [t2, y2] = ode15s(ode, tAdd:0.1:tEnd, y1(end,:),  opts);
    tAll = [t1; t2(2:end)];
    yAll = [y1; y2(2:end,:)];
end