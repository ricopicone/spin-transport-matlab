%---physical rules
global T ii g G D B_r c T1 rmax u0 T12 T13 mu_p mu_e B0 B_d kB temp freq duty rho_1_langevin rho_2_langevin rho_3_langevin  d_rho_1_langevin d_rho_2_langevin d_rho_3_langevin

% precise constants to be used
digits(50)
precise_constant = vpa(...
    -1*... % I think this is because mag grad is negative
    tanh(...
        vpa( ...
            1 + vpa(gamma_e/gamma_p) * vpa(Delta_3/Delta_2) ...
        )/vpa( vpa(gamma_e/gamma_p) * ( 1 + vpa(Delta_3/Delta_2) ) ) * ...
        vpa(...
            vpa(mu_e * B_d) / ...
            vpa( kB * temp ) ...
        ) ...
    ) ...
);

mu_p_kB = vpa(mu_p/kB);
mu_e_kB = vpa(mu_e/kB);

% Langevin (equilibrium solution)
rho_1_langevin= @(r) double(precise_constant);    % langevin nuclear polarization
rho_2_langevin= @(r)  tanh((B0-B_d*r)*mu_p_kB/temp);    % langevin nuclear polarization
rho_3_langevin= @(r)  tanh((B0-B_d*r)*mu_e_kB/temp);    % langevin electron polarization

% gradient of equilibrium solution
d_rho_1_langevin= @(r) 0;
d_rho_2_langevin= @(r)  double(-B_d*mu_p_kB/temp*(sech((B0-B_d*r)*mu_p_kB/temp))^2);
d_rho_3_langevin= @(r)  double(-B_d*mu_e_kB/temp*(sech((B0-B_d*r)*mu_e_kB/temp))^2);
