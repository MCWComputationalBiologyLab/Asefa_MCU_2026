function [beta_BPx_Ca, BPCa] = Ca_Buffering_in(Ca_x)

% Inputs:  Ca_x (M)
% Outputs: beta_BP_Ca, BPCa (M)

    K_CaBP = 2.5e-6;        % M
    BP_T   = 25e-3;         % M

    BPCa = (Ca_x .* BP_T) ./ (K_CaBP + Ca_x);
    beta_BPx_Ca = 1 + ((BP_T .* K_CaBP) ./ (K_CaBP + Ca_x).^2);

end
