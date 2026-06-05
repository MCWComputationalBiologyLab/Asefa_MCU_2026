function [beta_BP_Ca] = Ca_Buffering_in(Ca_x)
% K values corresponding to the dissociation constants of the associated reactions
K_CaBP1  = 1*2e-6;               % CaBP1 dissociation constants
K_CaBP2  = 1*1.7e-6;             % Ca6BP2 dissociation constants

% Total reactant concentration
BP1_T  = .010*15e-3;    % Total matrix Binding protein1 concentaration
BP2_T  = .010*15e-3;    % Total matrix Binding protein2 concentaration


% Calculating betas 
beta_BP_Ca1 = 1 + ((BP1_T*K_CaBP1)/(K_CaBP1+Ca_x)^2);
beta_BP_Ca2 = (36*(Ca_x^5)*BP2_T*K_CaBP2^6)/((K_CaBP2^6)+(Ca_x^6))^2;
beta_BP_Ca = 10*beta_BP_Ca1 + 10*beta_BP_Ca2;

end
