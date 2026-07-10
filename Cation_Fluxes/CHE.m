function flux = CHE(Ca_x, Ca_e, H_x, H_e)
% Ca-H exchanger flux (dimensionless scale)

% Parameters
K_Ca = 7.94e-06;   % M
K_H  = 2.83e-07;   % M

% Regulators & building blocks
Cx2 = Ca_x.^2;
He2 = H_e.^2;
Hx2 = H_x.^2;

regulator = Cx2 ./ (K_Ca^2 + Cx2);

% Mass-action term
numerator  = (Ca_x .* He2) - (Ca_e .* Hx2);

% Denominator (thermodynamic/kinetic form)
denominator = (K_H^2)*K_Ca + (K_H^2).*(Ca_x + Ca_e) + K_Ca.*(Hx2 + He2) + (Ca_x .* He2) + (Ca_e .* Hx2);

flux = regulator .* (numerator ./ denominator);
end
