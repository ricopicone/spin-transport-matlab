global T freq_p duty_p freq_e duty_e T_s_e T_s_p B1max_p B1max_e

%---simulation paramters
% Tsec =      T13sec*100;                  % sec ... max simulation time
Tsec = 10e-3;
T = tPFunc(Tsec); % dimensionless ... normalized max simulation time
rmaxnm = 100; % nm
rmax=rPFunc(rmaxnm*10^-9);        % dimensionless "rbar"
nr=         300;                % number of r spatial positions 
nt=         1000;               % number of time incrments

n_p = 1; % number of proton pulses
n_e = 1; % number of electron pulses

duty_p = .5; % duty cycle of proton pulses
duty_e = .5; % duty cycle of electron pulses

on_B1_p = 0; % nuclear B1 on?
on_B1_e = 1; % electron B1 on?

B1max_p = on_B1_p*B1max_p_nom;
B1max_e = on_B1_e*B1max_e_nom;

T_p = T/n_p; % dimensionless ... period of proton pulses
freq_p = 1/T_p; % dimensionless ... frequency of proton pulses
T_s_p_nice = 1/(4*freq_p); % nice period of sinusoidal transition for B1 pulse
T_s_p = min([... % period of sinusoidal transition for B1 pulse
  T_s_p_nice,...
  2*duty_p/freq_p,...
  2*(1-duty_p)/freq_p...
]); 
T_e = T/n_e; % dimensionless ... period of electron pulses
freq_e = 1/T_e; % dimensionless ... frequency of electron pulses
T_s_e_nice = 1/(4*freq_e); % nice period of sinusoidal transition for B1 pulse
T_s_e = min([... % period of sinusoidal transition for B1 pulse
  T_s_e_nice,...
  2*duty_e/freq_e,...
  2*(1-duty_e)/freq_e...
]); 


%---window
range = 1*rmax/2;                  % position plot range
edge = -1*rmax/2;                  % position plot offset

%---pulse parameters

fprintf('\n');
fprintf('Simulation Parameters\n');
fprintf('T                     = %g\n',T);
fprintf('nr                    = %g\n',nr);
fprintf('nt                    = %g\n',nt);
fprintf('rmax                  = %g\n',rmax);
fprintf('\n');