global T

%---simulation paramters
% Tsec =      T13sec*100;                  % sec ... max simulation time
Tsec = 10e-3;
T =         tPFunc(Tsec); % dimensionless ... normalized max simulation time
rmaxnm = 100; % nm
rmax=rPFunc(rmaxnm*10^-9);        % dimensionless "rbar"
nr=         300;                % number of r spatial positions 
nt=         1000;               % number of time incrments
freqHz = 1/(3*6*T13sec);                    % Hz ... frequency of electron pulses
freq = freqHz/tPFunc(1);        % dimensionless ... frequency of electron pulses
duty = 1;                     % duty cycle of electron pulses

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