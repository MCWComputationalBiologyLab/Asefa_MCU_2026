function [beta_BPe_Ca] = Ca_Buffering_out(Ca_e)

% Inputs:  Ca_e (M)
% Outputs: beta_BPe_Ca

    K_CaBP = 5e-6;   % M
    BP_T   = 1e-6;   % M
    beta_BPe_Ca = 1 + ((BP_T .* K_CaBP) ./ (K_CaBP + Ca_e).^2);

end
