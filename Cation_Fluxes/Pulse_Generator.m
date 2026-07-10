function J_ion_add = Pulse_Generator(t, Ca_add, t_add, dt_add)

t_rel = t - t_add;

% Gaussian pulse
J_ion_add = Ca_add * (1/(dt_add*sqrt(2*pi))) .* exp(-(t_rel).^2 ./ (2*dt_add^2));
end