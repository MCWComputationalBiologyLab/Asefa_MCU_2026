function flux = KHE(K_x, K_e, H_x, H_e)
% K+/H+ exchanger flux (dimensionless scale)

% Parameters
K_K = 1.98e-02;   % M
K_H = 7.29e-09;   % M

numerator   = (K_x .* H_e) - (K_e .* H_x);
denominator = K_K*K_H + K_H .* (K_x + K_e) + K_K .* (H_e + H_x) + (K_x .* H_e) + (K_e .* H_x);

flux = numerator ./ denominator;
end
