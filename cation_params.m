function p = cation_params()

    %% Initial conditions (13-state)
    % y = [K_e,K_x,Na_e,Na_x,Ca_e,Ca_x, x0,x1, z0,z1, S_O,S_I,S_R]
    p.IC0 = [120e-3, 150e-3, 10e-3, 10e-3, 100e-9, 50e-9, 1,0, 1,0, 1,0,0];

    %%  Compartment
    p.V_e = 500;        % cytosolic volume Ratio
    p.V_x = 1;          % matrix volume Ratio   
    p.pH_e = 7.2;
    p.pH_x = 7.4;
    p.Mg_e = 0.5e-3;    % external Mg2+ (M)
    p.dPsi = 180e-3;    % membrane potential (V)

    %% --- Transporter Vmax
    p.Vmax_Kleak = 1e-9;
    p.Vmax_KHE = 1e-9;
    p.Vmax_NHE = 12e-9;
    p.Vmax_NCE = 4.2e-9;
    p.Vmax_CHE = 1e-9;
    p.Vmax_MCU = 0.2e-6;

    %% MICU 3-state (WT) gating kinetics
    p.k1_on  = 0.1e12;  % on-rate for x0 -> x1   (M^-2 s^-1)
    p.k1_off = 1.0;     % off-rate x1 -> x0      (s^-1)
    p.k2_on  = 0.2e12;  % on-rate x1 -> x12      (M^-2 s^-1)
    p.k2_off = 1.0;     % off-rate x12 -> x1     (s^-1)

    %% MICU 3-state (WT) conductance weights ---
    p.g0  = 0.0;
    p.g1  = 0.5;
    p.g12 = 1.0;

    %% MICU2-KO 2-state gating (only used when mode == 'MICU2KO')
    p.kz_off  = 1.5;
    p.kz_on   = 1.5e12;  
    p.g0_2KO  = 0.0;
    p.g1_2KO  = 1.0;

    %% Constitutively-open gate value (MICU1KO / FORCED_OPEN)
    p.G_open = 1;

    %% Biphasic inactivation (S_O/S_I/S_R), Ca_x-dependent
    p.k_OI_max = 10;        % s^-1
    p.k_IO     = 1;         % s^-1
    p.k_IR_max = 10;        % s^-1
    p.k_RI     = 1;         % s^-1
    p.K_i      = 600e-9;    %
    p.K_a      = 730e-9;    % M
    p.n_hill   = 3;

    %%  Gating selector & inactivation setting
    p.mode = 'WT';          % 'WT' | 'MICU2KO' | 'MICU1KO' | 'FORCED_OPEN' | 'NO_MCU'
    p.use_inactivation = 0; % 0 -> P_open_x = 1 ; 1 -> P_open_x = S_O + S_R

    %%  Pulse settings 
    p.t_Ca_add = 700;
    p.Ca_Pulse = 0;
    p.pw_Ca    = 2.5;
    p.t_Na_add = 1e10;   % effectively "no Na pulse"
    p.Na_Pulse = 0;
    p.pw_Na    = 2.5;

end
