% T12Now = T12;
T12Now = inf;

p_init = PiconeTwoSpeciesTransport_ic(r);
% p_pulse3 =  rho_3_langevin(r)* 1/2 * (1 + tanh(1e1*(r + .9*rmax)));
% p_pulse3 =  rho_3_langevin(r)* 1 * (1-2*exp(-r.^2*1e-1));
% pulse_now = pulse(t,freq,duty);
% if pulse_now == 1 
%     T13Term = ( p(3) - p_pulse3 )/1e-18;
% else
%     T13Term = 0;
% end

% w = 4*rmax/20;
% r0 = -rmax + w/5;
% slope = 3e1;%abs(2*rho_3_langevin(0))/(w/20);
% r0 = 0;

% if ( ( r0 - w/2 ) < r < ( r0 + w/2 ) )
%     T13Now = T13/100;
% %     T13Term = ( p(3) - rho_3_langevin(r) * hole(r,r0,w,1e-1) )/(T13/10);
% else
% %     T13Term = ( p(3) - rho_3_langevin(r) )/T13;
%     T13Now = T13;
% end
T13Now = T13/2e4;

% T13Term = 0;
T13Term = ( p(3) - p_init(3) )/T13Now;
% T13Term = ( p(3) - rho_3_langevin(r) * hole(r,r0-1*w,w,2e-1) * hole(r,r0+1*w,w,2e-1) )/T13Now;

% T12Term = (p(2)-rho_2_langevin(r))/T12Now;
T12Term = 0;