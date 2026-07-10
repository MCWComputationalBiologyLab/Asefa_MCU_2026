function flux = Kleak(K_e, K_x, dPsi)

% Thermodynamic constants
F = 96484.6;       % Faraday's constant (J * V^-1 * mol^-1)
R = 8314e-3;       % Universal gas constant (J * mol^-1 * K^-1)
T = 298;           % Temperature (K)
FRT = F/(R*T);     % V^-1

% Kleak kinetic parameters
K_K        = 1.0e-3;   % M
beta_Kleak = 0.5;      % dimensionless

exp_e = exp(beta_Kleak .* FRT .* dPsi);          % Forward exponential term
exp_x = exp(-(1-beta_Kleak) .* FRT .* dPsi);     % Reverse exponential term

% K+ leak flux
flux = (1/K_K)*(K_e.* exp_e  - K_x .* exp_x);

end