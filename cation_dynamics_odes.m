function dydt = cation_dynamics_odes(t, y, p)
% CATION_DYNAMICS_ODES mitochondrial cation-dynamics ODE.

%   State vector (13 states):
%     y = [ K_e, K_x, Na_e, Na_x, Ca_e, Ca_x, x0, x1, z0, z1, S_O, S_I, S_R ]
%          1    2    3     4     5     6     7   8   9   10  11   12   13
%
%     x0,x1        : MICU 3-state (WT) gating occupancies (x12 = 1-x0-x1)
%     z0,z1        : MICU2-KO 2-state gating occupancies  (z1  = 1-z0)
%     S_O,S_I,S_R  : biphasic inactivation (open / inactivated / recovered)
%
%   Gating behaviour is selected by p.mode:
%     'WT'          : 3-state MICU gating (x-states)         x inactivation
%     'MICU2KO'     : 2-state gating (z-states)              x inactivation
%     'MICU1KO'     : constitutively open, G = p.G_open      x inactivation
%     'FORCED_OPEN' : constitutively open, G = p.G_open      (NO inactivation)
%     'NO_MCU'      : J_MCU = 0

    % pH -> [H+]
    H_e = 10^(-p.pH_e);
    H_x = 10^(-p.pH_x);

    % Pulses
    J_Ca_pulse = Pulse_Generator(t, p.Ca_Pulse, p.t_Ca_add, p.pw_Ca);
    J_Na_pulse = Pulse_Generator(t, p.Na_Pulse, p.t_Na_add, p.pw_Na);

    % States
    K_e = y(1);  K_x  = y(2);
    Na_e = y(3); Na_x = y(4);
    Ca_e = y(5); Ca_x = y(6);
    x0  = y(7);  x1   = y(8);  x2 = max(0, 1 - x0 - x1);
    z0  = y(9);  z1   = y(10); 
    S_O = y(11); S_I  = y(12);  S_R = y(13);

    %  Buffering factors 
    beta_BPe_Ca = Ca_Buffering_out(Ca_e);   % external buffer
    beta_BPx_Ca = Ca_Buffering_in(Ca_x);    % matrix buffer (scalar beta)

    %  Transporter fluxes 
    J_NCE = p.Vmax_NCE * NCE(Na_x, Na_e, Ca_x, Ca_e, p.dPsi);
    J_NHE = p.Vmax_NHE * NHE(Na_x, Na_e, H_x, H_e);
    J_CHE = p.Vmax_CHE * CHE(Ca_x, Ca_e, H_x, H_e);
    J_KHE = p.Vmax_KHE * KHE(K_x,  K_e,  H_x, H_e);
    J_Kleak = p.Vmax_Kleak * Kleak(K_x,  K_e, p.dPsi);

    %  Biphasic inactivation dynamics (Ca_x-dependent)
    k_OI = p.k_OI_max * Ca_x / (p.K_i + Ca_x);
    k_IR = p.k_IR_max * (Ca_x^p.n_hill) / (p.K_a^p.n_hill + Ca_x^p.n_hill);
    dS_O = -k_OI*S_O + p.k_IO*S_I;
    dS_I =  k_OI*S_O - (p.k_IO + k_IR)*S_I + p.k_RI*S_R;
    dS_R =  k_IR*S_I - p.k_RI*S_R;

    if p.use_inactivation
        P_open_x = S_O + S_R;
    else
        P_open_x = 1;
    end

    %  MICU gating
    C2  = max(Ca_e, eps)^2;
    dx0 = 0; dx1 = 0; dz0 = 0; dz1 = 0;

    switch upper(p.mode)
        case 'WT'
            dx0 = -p.k1_on*C2*x0 + p.k1_off*x1;
            dx1 =  p.k1_on*C2*x0 - (p.k1_off + p.k2_on*C2)*x1 + p.k2_off*x2;
            G_micu = min(1, max(0, p.g0*x0 + p.g1*x1 + p.g12*x2));
            G_gate = G_micu * P_open_x;

        case 'MICU2KO'
            dz0 = -p.kz_on*C2*z0 + p.kz_off*(1 - z0);
            dz1 = -dz0;
            G_micu = min(1, max(0, p.g0_2KO*z0 + p.g1_2KO*(1 - z0)));
            G_gate = G_micu * P_open_x;

        case 'MICU1KO'
            G_gate = p.G_open * P_open_x;          % inactivation applied

        case 'FORCED_OPEN'
            G_gate = p.G_open;                     % flat, no inactivation

        case 'NO_MCU'
            G_gate = 0;

    end
    G_gate = min(1, max(0, G_gate));

    %  MCU flux
    if strcmpi(p.mode, 'NO_MCU')
        J_MCU = 0;
    else
        J_MCU = p.Vmax_MCU * G_gate * MCU(Ca_x, Ca_e, p.Mg_e, p.dPsi);
    end

    % Mass balances
    dK_e  = (J_KHE - J_Kleak) / p.V_e;
    dK_x  = (J_Kleak - J_KHE) / p.V_x;
    dNa_e = (J_NHE - 3*J_NCE) / p.V_e + J_Na_pulse;
    dNa_x = (3*J_NCE - J_NHE) / p.V_x;  
    dCa_e = (1/beta_BPe_Ca) * (-2*J_MCU + J_NCE + J_CHE) / p.V_e + J_Ca_pulse/beta_BPe_Ca;
    dCa_x = (1/beta_BPx_Ca) * (2*J_MCU - J_NCE - J_CHE) / p.V_x;

    dydt = [ dK_e; dK_x; dNa_e; dNa_x; dCa_e; dCa_x; ...
             dx0;  dx1;  dz0;   dz1;   dS_O;  dS_I;  dS_R ];
end
