function flux = KHE(K_x, K_e, H_x, H_e)
% Vectorized K+/H+ exchanger flux (dimensionless scale)
% Inputs can be scalars or same-size arrays.

% Parameters
K_K = 1.9793e-02;   % M
K_H = 7.2865e-09;   % M

numerator   = (K_x .* H_e) - (K_e .* H_x);
denominator = K_K*K_H ...
            + K_H .* (K_x + K_e) ...
            + K_K .* (H_e + H_x) ...
            + (K_x .* H_e) + (K_e .* H_x);

denominator = denominator + eps(denominator);

flux = numerator ./ denominator;
end
