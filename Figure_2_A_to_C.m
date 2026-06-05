%% ------------------------------------------------------------------------- 
% MICU WT: Time-course gating at three [Ca2+] conditions with step at t_step
% States: x0 (no Ca bound), x1 (MICU1·2Ca), x2 = 1 - x0 - x1 (MICU2·4Ca)
% Output (Tot) = g1*x1 + g2*x2
% -------------------------------------------------------------------------

close all; clear; clc;

%  Simulation times 
t_end   = 50;           % Total simulation time (s)
t_step  = 0.15;         % Time when Ca2+ step is applied (seconds)

%  Ca2+ conditions 
Ca = [500e-9, 5e-6, 10e-6];

%  Gating Weights 
g1  = 0.5; 
g2  = 1;   

%  Figure 1 y-axis limit
ylm = [0 1];

%  Figure layout 
f  = figure('Position', [100, 200, 1150, 300], 'Color', 'w');   
tlo = tiledlayout(f, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

%  Loop panels 
for i = 1:numel(Ca)

    Ca_final = Ca(i);
    y0       = [1; 0];
    tspan    = [0, t_end];        

    % Solve ODEs
    [t, y] = ode15s(@(tt, yy) MICU12_WT(tt, yy, t_step, Ca_final), tspan, y0);

    % Extract states and total conductance
    x0  = y(:,1);
    x1  = y(:,2);
    x2  = 1 - x0 - x1;
    Tot = g1*x1 + g2*x2;

    % Plot
    nexttile(tlo); hold on; box on;

    h0 = plot(t, x0, '-', 'LineWidth', 1.5, 'Color', 'r');
    h1 = plot(t, x1, '-', 'LineWidth', 1.5, 'Color', 'g');
    h2 = plot(t, x2, '-', 'LineWidth', 1.5, 'Color', 'b');
    % h3 = plot(t, Tot, '-', 'LineWidth', 2.0, 'Color', 'k');

    if Ca_final < 1e-6
        ttl = sprintf('%.0f nM', Ca_final*1e9);
    else
        ttl = sprintf('%.1f \\muM', Ca_final*1e6);
    end

    title(ttl, 'FontWeight', 'normal', 'FontSize', 16);
    xlabel('Time (s)', 'FontSize', 16);

    if i == 1
        ylabel({'Fractional States'}, 'FontSize', 16);

        legend([h0 h1 h2], {'x_0', 'x_1', 'x_2'}, ...
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


%% Function: MICU_WT
function dydt = MICU12_WT(t, y, t_step, Ca)

    k1_off = 1;
    k2_off = 1;
    k1_on  = 1.1*k1_off / (3.5e-6)^2;
    k2_on  = 20*k2_off / (9.5e-6)^2;

    x0 = y(1);
    x1 = y(2);
    x2 = 1 - x0 - x1;

    if t < t_step
        C = 0;
    else
        C = Ca;
    end

    dx0_dt = -k1_on * (C^2) * x0 + k1_off * x1;
    dx1_dt =  k1_on * (C^2) * x0 ...
            - (k1_off + k2_on * (C^2)) * x1 ...
            + k2_off * x2;

    dydt = [dx0_dt; dx1_dt];
end