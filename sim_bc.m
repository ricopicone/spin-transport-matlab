% --------------------------------------------------------------
function [pl,ql,pr,qr] = sim_bc(xl,ul,xr,ur,t)
global u0 c g G

%---pinned
% pl = [ul(1); ul(2)-u0(2); ul(3)-u0(3)];                               
% ql = [0;0;0];                                  
% pr = [ur(1); ur(2)-u0(2); ur(3)-u0(3)];                         
% qr = [0;0;0];      

% %---zero gradient
% pl = [0; 0; 0];                               
% ql = [1;1;1];                                  
% pr = [0; 0; 0];                         
% qr = [1;1;1];      

%---zero j
pl = [0; c*(1-ul(2)^2)*ul(1); -c*g*(1-ul(3)^2)*ul(1)];                               
ql = [1/(1+G);1;1/G];                                  
pr = [0; c*(1-ur(2)^2)*ur(1); -c*g*(1-ur(3)^2)*ur(1)];                         
qr = [1/(1+G);1;1/G];  

end