%---physical rules
global T ii g G D B_r c T1 rmax u0 T12 T13 mu_p mu_e B0 B_d kB temp freq duty

% Langevin
rho_1_langevin= @(r)  ( 1 + gamma_e/gamma_p * Delta_3/Delta_2 ) * mu_p * B_d / ( kB * temp );    % langevin nuclear polarization
rho_2_langevin= @(r)  tanh((B0-B_d*r)*mu_p/kB/temp);    % langevin nuclear polarization
rho_3_langevin= @(r)  tanh((B0-B_d*r)*mu_e/kB/temp);    % langevin electron polarization