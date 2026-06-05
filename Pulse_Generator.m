function J_ion_add = Pulse_Generator(t, Ca_add, t_add, dt_add)

% % If before addition time, return 0
% if t < t_add || (t > t_add + 10) % Cut off after 10 seconds
%     J_ion_add = 0;
%     return;
% end

% Calculate time relative to addition point
% t_rel = t - t_add;
% 
%     % Gaussian pulse
%     J_ion_add = Ca_add*(1/(dt_add*sqrt(2*pi)))*exp(-(t_rel)^2/(2*dt_add^2));
% end

t_rel = t - t_add;

% Gaussian pulse
J_ion_add = Ca_add * (1/(dt_add*sqrt(2*pi))) .* ...
            exp(-(t_rel).^2 ./ (2*dt_add^2));
end