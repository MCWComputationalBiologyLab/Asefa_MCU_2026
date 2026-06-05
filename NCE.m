function J_NCE = NCE(Na_x, Na_e, Ca_x, Ca_e)

%% Standard parameter values
F    = 96484.6;    % Faraday's constant (C/mol)
R    = 8.314;      % Universal gas constant (J/mol/K)
T    = 298;        % Temperature (K)
RT   = R*T;        % J/mol

% Optimal parameters
K_Na = 8.2e-3;         % M (Paucek and Jaburek, 2004)
K_Ca = 2.1e-6;         % M (Paucek and Jaburek, 2004)
V_NCE = 0.024*0.176;   % umol/mg/s (estimated to fit the data)

dPsi = 190e-3;               % V

n_Ca  = 0.5;                  % fractional coefficient of membrane potential

%% NaCa Exchanger flux with Ca_e addition and Na_x exclusion (umol/mg/s)
% For Ca to go in and Na to come out, V_NCE is negative.
J_NCE = -V_NCE * (exp(-n_Ca*F*dPsi/RT) * (Na_e.^3 .* Ca_x)/(K_Na^3*K_Ca) - ...
         exp((1-n_Ca)*F*dPsi/RT) * (Na_x.^3 .* Ca_e)/(K_Na^3*K_Ca)) ./ ...
        (1 + Na_e.^3/K_Na^3 + Ca_x/K_Ca + (Na_e.^3 .* Ca_x)/(K_Na^3*K_Ca) + ...
           Na_x.^3/K_Na^3 + Ca_e/K_Ca + (Na_x.^3 .* Ca_e)/(K_Na^3*K_Ca));
end
