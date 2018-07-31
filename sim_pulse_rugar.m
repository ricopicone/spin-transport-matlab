% T12Now = T12;
T12Now = inf;

p_init = PiconeTwoSpeciesTransport_ic(r);
p_pulse3 = p_init(3);
pulse_now = pulse(t,freq,duty);
% if pulse_now == 1 
%     T13Term = ( p(3) - p_pulse3 )/T13;
% else
%     T13Term = ( p(3) - rho_3_langevin(r) )/T13;
% end
% T13Term = ( p(3) - p_init(3) )/T13;
% T13Term = 0;
w = 2*rmax/20;
vwave = -(2*(rmax - w)/T)/4;
% r0 = vwave*t - rmax/2;
% r0 = -rmax + w;
r0 = 0;
% if ( ( r0 - w/2 ) < r < ( r0 ) )
%     T13Now = T13/20;
% %     T13Term = ( p(3) - rho_3_langevin(r) * hole(r,r0,w,1e-1) )/(T13/10);
% else
% %     T13Term = ( p(3) - rho_3_langevin(r) )/T13;
%     T13Now = T13;
% end
T13Now = T13;

% T13Term = ( p(3) - rho_3_langevin(r) * hole(r,r0,w,2e-1) )/T13Now;
T13Term = ( p(3) - rho_3_langevin(r) * hole(r,r0-1*w,w,2e-1) * hole(r,r0+1*w,w,2e-1) )/T13Now;

% T12Term = (p(2)-rho_2_langevin(r))/T12Now;
T12Term = 0;