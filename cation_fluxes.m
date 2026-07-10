function flux = cation_fluxes(t, Y, p)
% CATION_FLUXES  Recompute fluxes and the effective gate from a trajectory.
%
%   flux = cation_fluxes(t, Y, p)
%
%   Given a simulated trajectory (t, Y) produced by cation_dynamics_odes with
%   parameter struct p, this recomputes all transporter fluxes and the effective
%   MCU open gate using the SAME p. Because it reuses p (including p.Vmax_MCU and
%   p.mode), the plotted J_MCU is guaranteed to match the dynamical J_MCU that
%   drove the simulation.
%
%   Returns struct flux with fields:
%     J_MKU, J_KHE, J_NHE, J_NCE, J_CHE, J_MCU 
%     G        effective open gate that multiplies MCU
%     G_micu   MICU conductance before inactivation  

    H_e = 10^(-p.pH_e);
    H_x = 10^(-p.pH_x);

    K_e = Y(:,1);  K_x  = Y(:,2);
    Na_e= Y(:,3);  Na_x = Y(:,4);
    Ca_e= max(Y(:,5), 0);
    Ca_x= max(Y(:,6), 0);
    x0  = Y(:,7);  x1   = Y(:,8);  x12 = max(0, 1 - x0 - x1);
    z0  = Y(:,9);
    S_O = Y(:,11); S_R  = Y(:,13);

    if p.use_inactivation
        P_open_x = S_O + S_R;
    else
        P_open_x = ones(size(S_O));
    end

    switch upper(p.mode)
        case 'WT'
            G_micu = min(1, max(0, p.g0.*x0 + p.g1.*x1 + p.g12.*x12));
            G = G_micu .* P_open_x;
        case 'MICU2KO'
            z1 = max(0, 1 - z0);
            G_micu = min(1, max(0, p.g0_2KO.*z0 + p.g1_2KO.*z1));
            G = G_micu .* P_open_x;
        case 'MICU1KO'
            G_micu = p.G_open .* ones(size(x0));
            G = G_micu .* P_open_x;
        case 'FORCED_OPEN'
            G_micu = p.G_open .* ones(size(x0));
            G = G_micu;                          % flat, no inactivation
        case 'NO_MCU'
            G_micu = zeros(size(x0));
            G = G_micu;
        otherwise
            error('cation_fluxes:badMode', 'Unknown p.mode = "%s".', p.mode);
    end
    G = min(1, max(0, G));
    n = numel(t);
    J_MCU = zeros(n,1);  J_NCE = zeros(n,1);  J_NHE = zeros(n,1); 
    J_CHE = zeros(n,1);  J_KHE = zeros(n,1);  J_Kleak = zeros(n,1); 

    for k = 1:n
        J_NCE(k) = p.Vmax_NCE * NCE(Na_x(k), Na_e(k), Ca_x(k), Ca_e(k), p.dPsi);
        J_NHE(k) = p.Vmax_NHE * NHE(Na_x(k), Na_e(k), H_x, H_e);
        J_CHE(k) = p.Vmax_CHE * CHE(Ca_x(k), Ca_e(k), H_x, H_e);
        J_KHE(k) = p.Vmax_KHE * KHE(K_x(k),  K_e(k),  H_x, H_e);
        J_Kleak(k) = p.Vmax_Kleak * Kleak(K_x(k),  K_e(k), p.dPsi);
        if strcmpi(p.mode, 'NO_MCU')
            J_MCU(k) = 0;
        else
            J_MCU(k) = p.Vmax_MCU * G(k) * MCU(Ca_x(k), Ca_e(k), p.Mg_e, p.dPsi);
        end
    end

    flux.J_MCU = J_MCU;  flux.J_NCE = J_NCE;    flux.J_NHE = J_NHE;  
    flux.J_CHE = J_CHE;  flux.J_KHE = J_KHE;    flux.J_Kleak = J_Kleak;  
    flux.G = G;          flux.G_micu = G_micu;  flux.P_open_x = P_open_x;
end
