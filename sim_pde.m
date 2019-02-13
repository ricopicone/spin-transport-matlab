% --------------------------------------------------------------
function [cm,f,s] = sim_pde(r,t,p,DpDx)

% constants and parameters
global ii g G D c rho_2_langevin rho_3_langevin T1 mu_p mu_e B0 Bd_2 Bd_3 kB temp T12 T13 T22 T23 freq duty rmax sim T B1max_p B1max_e freq_e duty_e T_s_e freq_p duty_p T_s_p gamma_p gamma_e

cm = [1; 1; 1];                                  
f = [...
    1+G; ...
    1; ...
    G ...
].*DpDx;

% separation terms
% s = [...
%     -c^2/(1+D)*( (1-p(2)^2) + G*D*g^2*(1-p(3)^2) )*atanh(p(1))-c/(1+D)*(DpDx(2)-G*D*g*DpDx(3)); ...
%     c*( (1-p(2)^2)/(1-p(1)^2)*DpDx(1) - 2*atanh(p(1))*p(2)*DpDx(2) ); ...
%     -c*G*g*( (1-p(3)^2)/(1-p(1)^2)*DpDx(1) - 2*atanh(p(1))*p(3)*DpDx(3) ) ...
% ];  
s = [0;0;0]; % no cross terms

% Bloch dynamics sans relaxation
% - amplitude modulation
% w1_2 = -gamma_p*B1(t,freq_p,duty_p,T_s_p,B1max_p);
% w1_3 = -gamma_e*B1(t,freq_e,duty_e,T_s_e,B1max_e);
% - frequency modulation
% % delta_tilde_2 = gamma_p*Bd_2*(r-rs);
% % delta_tilde_3 = gamma_e*Bd_3*(r-rs);
% delta_tilde_2 = gamma_p*Bd_2*r;
% delta_tilde_3 = gamma_e*Bd_3*r;
% - modulation enters through the taus
% tau_2 = 1/(1/T12+T22*w1_2^2/(1+T22^2*delta_tilde_2^2));
% tau_3 = 1/(1/T13+T23*w1_3^2/(1+T23^2*delta_tilde_3^2));
% s = s + [ ...
%     0; ...
%     -p(2)/tau_2; ...
%     -p(3)/tau_3 ...
% ];

% Bloch relaxation dynamics <--this slows the solver

% s = s + [ ...
%     0; ...
%     rho_2_langevin(r)/T12; ...
%     rho_3_langevin(r)/T13 ...
% ];

% print

ii = ii+1;
if mod(ii,10000)==0
    fprintf('%d  %g\n',ii,t)
end

end