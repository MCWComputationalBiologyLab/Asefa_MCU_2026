function flux = MKU(K_x, K_e)
% Vectorized monovalent cation uniporter–like flux (dimensionless scale)

% Parameters
X_MKU = 4.1456e-02;  % activity coefficient (arbitrary scale)
K_K   = 1.5699e-03;  % M

% Thermo
z  = 1;
F  = 96484.6;
R  = 8.314;
T  = 298;
FRT = F/(R*T);

% Voltage dependence
dPsi = 190e-3;                 % V (constant here)
nH   = 1;                      % Hill
dPhi = z.*dPsi.*FRT;
sinhx = sinh(dPhi./nH);
sinhx = sign(sinhx).*max(abs(sinhx), realmin);
Fact = ((dPhi./nH)./sinhx).^nH;

% Partitioning terms
den   = (K_K + K_x + K_e);
den   = den + eps(den);
Ge    = K_e ./ den;
Gx    = K_x ./ den;

flux = X_MKU .* Fact .* (Ge .* exp(dPhi) - Gx .* exp(-dPhi));
end
