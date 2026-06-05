function [beta_EGTA_Ca] = Ca_Buffering_out(Ca_e, pH_e)
% pK values corresponding to the dissociation constants of the associated reactions
pK_HEGTA   = 9.58;           % HEGTA(3-)  = H(+) + EGTA(4-)
pK_H2EGTA  = 8.96;           % H2EGTA(2-) = H(+) + HEGTA(3-)
pK_CaEGTA  = 1*10.97;          % CaEGTA(2-) = Ca(2+) + EGTA(4-)

% K values corresponding to the dissociation constants of the associated reactions
K_HEGTA  = 10^(-pK_HEGTA);     % HEGTA dissociation constants
K_H2EGTA = 10^(-pK_H2EGTA);    % H2EGTA dissociation constants
K_CaEGTA = 10^(-pK_CaEGTA);    % CaEGTA dissociation constants

% Total reactant concentration
EGTA_T = 2.5*40e-6;    % Total buffer EGTA concentaration
H_e = 10^-pH_e;    % free H+ external concentration

% Calculating binding polynomials for EGTA
P_EGTA = 1 + (H_e/K_HEGTA)+((H_e^2)/(K_HEGTA*K_H2EGTA))+(Ca_e/K_CaEGTA);

% Apparent dissociation constants (Kapp) for CaEGTA
K_app_CaEGTA = (K_CaEGTA*P_EGTA) - Ca_e;  

% Calculating beta 
beta_EGTA_Ca = 1 + ((EGTA_T*K_app_CaEGTA)/((K_app_CaEGTA+Ca_e)^2));
end
