% --------------------------------------------------------------
function [pl,ql,pr,qr] = sim_bc(xl,ul,xr,ur,t)
global u0 c g G rho_1_langevin rho_2_langevin rho_3_langevin d_rho_1_langevin d_rho_2_langevin d_rho_3_langevin

%---pinned
% pl = ul - [rho_1_langevin(xl);rho_2_langevin(xl);rho_3_langevin(xl)];
% ql = [0;0;0];
% pr = ur - [rho_1_langevin(xr);rho_2_langevin(xr);rho_3_langevin(xr)];
% qr = [0;0;0];      

% %---zero gradient
% pl = [0; 0; 0];                               
% ql = [1;1;1];                                  
% pr = [0; 0; 0];                         
% qr = [1;1;1];      

%---zero j
% pl = [0; -c*(1-ul(2)^2)*atanh(ul(1)); -c*G*g*(1-ul(3)^2)*atanh(ul(1))];                               
% ql = [-1;-1;-1];                                  
% pr = [0; -c*(1-ur(2)^2)*atanh(ur(1)); -c*G*g*(1-ur(3)^2)*atanh(ur(1))];                         
% qr = [-1;-1;-1];     

%---equil j                             
ql = [-1/(1+G);-1;-1/G];                                                         
qr = [-1/(1+G);-1;-1/G];  
pl = [d_rho_1_langevin(xl);d_rho_2_langevin(xl);d_rho_3_langevin(xl)];
pr = [d_rho_1_langevin(xr);d_rho_2_langevin(xr);d_rho_3_langevin(xr)];

end