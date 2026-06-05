%% -------------------------------------------------------------------------
% MICU2-KO (2-state model): Time-course gating at three Ca2+ levels
% This script simulates the time-dependent behavior of the MICU2-KO system,
% which lacks the MICU2 binding step. The remaining dynamics correspond to 
% a single MICU1-dependent binding transition.
% -------------------------------------------------------------------------

close all; clear; clc;

%  Simulation times 
t_end   = 50;           % Total simulation time (s)
t_step  = 0.15;         % Time when Ca2+ step is applied (seconds)

%  Ca2+ conditions 
Ca = [500e-9, 5e-6, 10e-6];

%  Gating weights 
g1 = 1;    % high open probability upon Ca2+ binding

%  Figure y-limit
ylm = [0 1];

%  Figure layout 
f  = figure('Position', [100, 200, 1150, 300], 'Color', 'w');   
tlo = tiledlayout(f, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

%  Loop panels 
for i = 1:numel(Ca)
    Ca_final = Ca(i); 
    y0       = [1; 0];              % initial state: all in x0 (resting)
    tspan    = [0, t_end];

    %  Solve ODEs 
    [t, y] = ode15s(@(tt, yy) MICU2_KO(tt, yy, t_step, Ca_final), tspan, y0);

    % Extract states and total conductance 
    x0  = y(:,1);      % unbound fraction
    x1  = y(:,2);      % bound fraction
    Tot = g1*x1;       % weighted sum (P_open,c)
 
    %  Plot
    nexttile(tlo); hold on; box on;
    h0 = plot(t, x0, '-', 'LineWidth', 1.5, 'Color', 'r');  % x0 (red)
    h1 = plot(t, x1, '-', 'LineWidth', 1.5, 'Color', 'g');  % x1 (green)
    % h2 = plot(t, Tot, '-', 'LineWidth', 2.0, 'Color', 'k');  % total gating output

    if Ca_final < 1e-6
        ttl = sprintf('%.0f nM', Ca_final*1e9);
    else
        ttl = sprintf('%.1f \\muM', Ca_final*1e6);
    end

    title(ttl, 'FontWeight', 'normal', 'FontSize', 16);
    xlabel('Time (s)', 'FontSize', 16);

    % Show ylabel and legend only on the first subplot
    if i == 1
        ylabel({'Fractional States'}, 'FontSize', 16);

        legend([h0 h1], {'x_0', 'x_1'}, ...
            'Box', 'off', ...
            'FontSize', 17, ...
            'Location', 'west');
    end

    xlim([0 2]);
    ylim(ylm);
    xticks(0:1:2);
    yticks(0:0.5:1);

    set(gca, 'FontSize', 16, 'LineWidth', 1.0);
end


%% Function: MICU2_KO
function dydt = MICU2_KO(t, y, t_step, Ca)

    %  Kinetic Parameters 
    k1_on  = 1.5 / (1.0e-6)^2;    % MICU1 binding rate (M^-2 s^-1)
    k1_off = 1.5;                 % MICU1 unbinding rate (s^-1)

    %  States 
    x0 = y(1);                    % no Ca2+ bound
    x1 = y(2);                    % MICU1-Ca2+ bound

    %  Ca2+ Step 
    if t < t_step
        C = 0;
    else
        C = Ca;
    end

    %  Differential Equations 
    dx0_dt = -k1_on * (C^2) * x0 + k1_off * x1;
    dx1_dt =  k1_on * (C^2) * x0 - k1_off * x1;

    dydt = [dx0_dt; dx1_dt];
end