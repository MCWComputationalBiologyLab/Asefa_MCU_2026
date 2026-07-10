function flux = MCU(Ca_x, Ca_e, Mg_e, dPsi)
% Thermodynamic constants
Z_Ca = 2;          % Valency of Ca (uniteless)
F = 96484.6;       % Faraday's constant (J * V^-1 * mol^-1) (Coul * mol^-1)
R = 8314e-3;       % Universal gas constant (J * mol^-1 * K^-1)
T = 298;           % Temperature (K)
FRT = F/(R*T);     % V^-1 

% Calcium uniporter model parameter values
K_Ca  = 4e-6;    % Binding constant for Ca on MCU (M)
K_Mg  = 0.3e-3;  % Binding constant for Mg on MCU (M)
alpha = 55;      % Scaling for Mg2+ inhibition type
nH = 2.7;        % Hill coefficient for membrane potential dependency

%  Eyring-type voltage term 
dPhi = Z_Ca*dPsi*FRT;
if (nH == 0 || dPhi == 0)
    E_dPhi = 1;
else
    E_dPhi = ((dPhi/nH)/sinh(dPhi/nH))^nH;  % Voltage-dependent correction
end

beta  = 0.5*(1 + log(E_dPhi)/dPhi);         % Symmetry factor (dimensionless)
exp_e = exp(+2*beta*dPhi);                  % Forward exponential term
exp_x = exp(-2*(1-beta)*dPhi);              % Reverse exponential term


K_Ce = K_Ca;   K_Cx = K_Ce * 10;     % Matrix side weaker binding

num   = ((Ca_e./K_Ce).^2 * exp_e - (Ca_x./K_Ce).^2 * exp_x);
denom = 1 + (Ca_e./K_Ce).^2 + Mg_e./K_Mg + (Ca_e./K_Ce).^2 .* (Mg_e./(alpha*K_Mg)) + (Ca_x./K_Cx).^2;

%  Final flux expression 
flux = (num ./ denom);  % nmol/mg/s

end
