function flux = NHE(Na_x, Na_e, H_x, H_e)
% Vectorized Na+/H+ exchanger flux (dimensionless scale)

K_Na = 2.1999e-02;   % M
K_H  = 9.9825e-08;   % M

% Regulation by matrix H+
% Guard (avoid huge factors if H_x ~ 0): add tiny epsilon in denom
regulation_factor = (H_x.^2) ./ max(K_H .* (H_x + K_H), realmin);

numerator   = (Na_x .* H_e) - (Na_e .* H_x);
denominator = K_Na .* (H_e + H_x) + (Na_x .* H_e) + (Na_e .* H_x);
denominator = denominator + eps(denominator);

flux = regulation_factor .* (numerator ./ denominator);
end
