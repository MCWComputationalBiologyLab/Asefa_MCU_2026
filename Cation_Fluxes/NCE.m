function J_NCE = NCE(Na_x, Na_e, Ca_x, Ca_e, dPsi)

%% Standard parameter values
F    = 96484.6;    % Faraday's constant (C/mol)
R    = 8.314;      % Universal gas constant (J/mol/K)
T    = 298;        % Temperature (K)
RT   = R*T;        % J/mol

% Optimal parameters
K_Na = 8.2e-3;         % M (Paucek and Jaburek, 2004)
K_Ca = 2.1e-6;         % M (Paucek and Jaburek, 2004)

%% NaCa Exchanger flux

numerator =  exp(0.5*F*dPsi/RT) .* (Na_e.^3 .* Ca_x) ./ (K_Na^3*K_Ca) ...
           - exp(-0.5*F*dPsi/RT) .* (Na_x.^3 .* Ca_e) ./ (K_Na^3*K_Ca);

denominator = (1 + Na_e.^3/K_Na^3 + Ca_x/K_Ca + (Na_e.^3 .* Ca_x)/(K_Na^3*K_Ca) + ...
           Na_x.^3/K_Na^3 + Ca_e/K_Ca + (Na_x.^3 .* Ca_e)/(K_Na^3*K_Ca));

J_NCE = numerator ./ denominator;

end