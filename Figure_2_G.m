%% -------------------------------------------------------------------------
% MICU gating: MICU12-WT (3-state) vs MICU2-KO (2-state) vs MICU1-KO (assumed G=1)
% Steady-state P_open,c vs cytosolic [Ca2+] from 100 nM to 100 µM
% WT states:  x0 (no Ca), x1 (MICU1-2Ca), x2 = 1 - x0 - x1 (MICU2-4Ca)
% KO2 states: x0 (no Ca), x1 (MICU1-Ca-bound)
% Conductance (P_open,c):
%   WT:    G = g1_WT*x1 + g2_WT*x2
%   KO2:   G = g1_KO*x1
%   KO1:   G = 1 (assumed constitutively open)
% Kinetics chosen via k_on = k_off / Kd^2 (two Ca2+ per step) and trial-tuned.
% -------------------------------------------------------------------------

close all; clear; clc;

Z = [0.5296557  0.0177384
0.6542395  0.0443459
0.7091955  0.0133038
0.7776771  0.0288248
0.7927561  0.0576497
0.8693066  0.0820399
0.9717322  0.0532151
0.9569172  0.0687361
0.9068210  0.0687361
1.1072879  0.0133038
1.2330087  0.0354767
1.1956985  0.1019956
1.3061317  0.0376940
1.3263532  0.0509978
1.4940632  0.0598670
1.5525641  0.0931264
1.4883357  0.0953437
1.2666089  0.1685144
1.6765274  0.0997783
1.6959570  0.1263858
1.1820001  0.0864745
1.7759532  0.3082040
2.0868436  0.2128603
2.7410813  0.3414634
2.4147715  0.4722838
3.1354791  0.6141907
3.4382485  0.6297118
3.4780950  0.6518847
3.2457601  0.6651885
3.6142818  0.6917960
4.1343191  0.6674058
4.1661999  0.6784922
8.1906447  0.8159645
9.1206055  0.8470067
16.6689628  0.9268293
18.3488984  0.9645233
22.8396386  0.9667406
27.7817834  0.9689579
26.8378431  0.9068736
26.5303781  0.8492239
29.5426213  0.9733925
23.9169562  0.9667406
24.7581631  0.9800443
18.0691519  0.9401330
21.8108478  0.9401330
22.8396386  0.9157428
24.6632532  0.9223947
2.0708745  0.3414634
2.0948742  0.3725055
2.6581377  0.3858093
2.7728482  0.4146341
1.4322550  0.1884701
1.2569164  0.1485588
1.1820001  0.0620843
0.9138138  0.0487805
7.2433633  0.9135255
7.8217045  0.9556541
8.4138445  0.8337029
8.5769873  0.8381375
14.6283266  0.9733925
14.4053038  0.9800443
17.3883031  0.9423503
16.1645688  0.9401330
29.5426213  0.9955654
31.1746669  0.9733925
32.3953286  0.9689579
];
% Gating weights
% MICU12-KO
g1_WT  = 0.45;   % partial activation (MICU1-bound)
g2_WT = 0.96;    % full activation (MICU2-bound)

% MICU2-KO
g1_KO  = 0.98;   % active when MICU1 is Ca-bound

%  [Ca2+] sweep (100 nM -> 100 µM) 
C_range = logspace(log10(100e-9), log10(100e-6), 100);

T_max = 5.0;          % integrate long enough to reach steady state
y0_WT = [1; 0];       % WT starts fully in x0
y0_KO = [1; 0];       % KO starts fully in x0


opts  = odeset('RelTol',1e-9,'AbsTol',1e-12,'NonNegative',[1 2]);

%  Preallocation 
G_WT   = zeros(size(C_range));
x0_WT  = zeros(size(C_range));  % MICU12-WT x0
x1_WT  = zeros(size(C_range));  % MICU12-WT x1
x2_WT  = zeros(size(C_range));  % MICU12-WT x12

G_KO   = zeros(size(C_range));
x0_KOa = zeros(size(C_range));  % MICU2-KO x0
x1_KOa = zeros(size(C_range));  % MICU2-KO x1

%  Parameter
for i = 1:numel(C_range)
    Ca = C_range(i);

    %  MICU12-WT (3-state) 
    [~, Ywt] = ode15s(@(t,y) MICU_WT(t, y, Ca), [0 T_max], y0_WT, opts);
    y_ss_wt = Ywt(end, :).';        
    y_ss_wt = max(y_ss_wt, 0);        
    x0 = y_ss_wt(1);
    x1 = y_ss_wt(2);
    x2 = 1 - x0 - x1;     

    x0_WT(i)  = x0;
    x1_WT(i)  = x1;
    x2_WT(i) = x2;
    G_WT(i)   = g1_WT*x1 + g2_WT*x2;

    %  MICU2-KO (2-state) 
    [~, Yko] = ode15s(@(t,y) MICU2_KO(t, y, Ca), [0 T_max], y0_KO, opts);
    y_ss_ko = Yko(end, :).';
    y_ss_ko = max(y_ss_ko, 0);
    x0k = y_ss_ko(1);
    x1k = y_ss_ko(2);

    x0_KOa(i) = x0k;
    x1_KOa(i) = x1k;
    G_KO(i)   = g1_KO*x1k;
end

%  Plotting 
figure; hold on; box on;
plot(C_range*1e6, G_WT, 'k-', 'LineWidth', 2);                 % WT (black)
plot(C_range*1e6, G_KO, 'b-', 'LineWidth', 2);                 % MICU2-KO (blue)
plot(C_range*1e6, ones(size(C_range)), 'r-', 'LineWidth', 2);  % MICU1-KO (assumed G=1)
scatter(Z(:,1), 1*Z(:,2), 30, 'k', 'filled', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
legend({'WT(MICU1/2)','MICU2-KO','MICU1-KO','WT Data'}, 'Location','southeast','Box','off');
xlabel('[Ca^{2+}]_c (\muM)', 'FontSize', 16);
ylabel('P_{open,c}',       'FontSize', 16);
set(gca, 'XScale', 'log', 'FontSize', 20, 'LineWidth', 1.5, ...
         'YLim', [0 1.05], 'XLim', [1e-1 1e2], 'YTick', 0:0.2:1);


% Function: MICU12-WT ODE (3-state x0 <-> x1 <-> x2)
function dydt = MICU_WT(~, y, Ca)
    % k1_off = 1.5;                 % MICU1 unbinding (s^-1)
    % k2_off = 15e-1;                 % MICU2 unbinding (s^-1)
    % k1_on  = 1.25*k1_off / (3.5e-6)^2; % MICU1 binding (M^-2 s^-1)
    % k2_on  = 10*k2_off / (9.5e-6)^2; % MICU2 binding (M^-2 s^-1)

    k1_off = 1;                 % MICU1 unbinding (s^-1)
    k2_off = 1;                 % MICU2 unbinding (s^-1)
    k1_on  = 0.1e12;%1.1*k1_off / (3.5e-6)^2; % MICU1 binding (M^-2 s^-1)
    k2_on  = 0.2e12;%20*k2_off / (9.5e-6)^2; % MICU2 binding (M^-2 s^-1)
    x0 = y(1);
    x1 = y(2);
    x2 = 1 - x0 - x1;   

    dx0_dt = -k1_on * (Ca^2) * x0 + k1_off * x1;
    dx1_dt =  k1_on * (Ca^2) * x0 - (k1_off + k2_on*(Ca^2)) * x1 + k2_off * x2;

    dydt = [dx0_dt; dx1_dt];
end

% Function: MICU2-KO ODE (2-state: x0 <-> x1)
function dydt = MICU2_KO(~, y, Ca)
    % k1_off = 1.5;                 % MICU1 unbinding (s^-1)
    % k1_on  = 1.5 / (1.0e-6)^2;    % MICU1 binding (M^-2 s^-1)
    k1_off = 1.5;                 % MICU1 unbinding (s^-1)
    k1_on  = 1.5 / (1.0e-6)^2;    % MICU1 binding (M^-2 s^-1)

    x0 = y(1);
    x1 = y(2);

    dx0_dt = -k1_on*(Ca^2)*x0 + k1_off*x1;
    dx1_dt =  k1_on*(Ca^2)*x0 - k1_off*x1;
    dydt = [dx0_dt; dx1_dt];
end
