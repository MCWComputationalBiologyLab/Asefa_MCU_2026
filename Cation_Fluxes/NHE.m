function flux = NHE(Na_x, Na_e, H_x, H_e)
% Na+/H+ exchanger flux

K_Na = 2.2e-02;   % M
K_H  = 1.0e-07;   % M

% Regulation by matrix H+

regulation_factor = (H_x.^2) ./ max(K_H .* (H_x + K_H));

numerator   = (Na_x .* H_e) - (Na_e .* H_x);
denominator = K_Na .* (H_e + H_x) + (Na_x .* H_e) + (Na_e .* H_x);

flux = regulation_factor .* (numerator ./ denominator);
end
