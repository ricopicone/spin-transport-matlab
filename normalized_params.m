%---normalized parameters
fprintf('\n');
fprintf('Normalized Parameters\n');
g = gamma_e/gamma_p;    
G = Gamma_3/Gamma_2;    
D = Delta_3/Delta_2;    
B_r = 1;    
c = B_r*(1+D)/(1+g*D);  
% T1 = inf;  

nd_t = @(t) Gamma_2*(grad/B_d)^2*t; % nondimensionalize time
nd_r = @(r) grad/B_d * r; % nondimensionalize space

% print
% fprintf('g =   gamma_e/gamma_p = %g\n',g);
% fprintf('G =   Gamma_3/Gamma_2 = %g\n',G);
% fprintf('D =   Delta_3/Delta_2 = %g\n',D);
% fprintf('B_r     = %g\n',B_r);
% fprintf('c = B_r*(1+D)/(1+g*D) = %g\n',c);
% fprintf('T1 fake = %g\n',T1);