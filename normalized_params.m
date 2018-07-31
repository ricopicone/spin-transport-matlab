%---normalized parameters
fprintf('\n');
fprintf('Normalized Parameters\n');
g=          gamma_e/gamma_p;    fprintf('g =   gamma_e/gamma_p = %g\n',g);
G=          Gamma_3/Gamma_2;    fprintf('G =   Gamma_3/Gamma_2 = %g\n',G);
D=          Delta_3/Delta_2;    fprintf('D =   Delta_3/Delta_2 = %g\n',D);
B_r=        1;                  fprintf('B_r                   = %g\n',B_r);
c=          B_r*(1+D)/(1+g*D);  fprintf('c = B_r*(1+D)/(1+g*D) = %g\n',c);
T1=         inf;                fprintf('T1 fake               = %g\n',T1);