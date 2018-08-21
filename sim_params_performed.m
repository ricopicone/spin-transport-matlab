%---simulation paramters
Tsec =      1e-2;                  % sec ... max simulation time
T =         tPFunc(Tsec); % dimensionless ... normalized max simulation time
rmaxnm = 100;                   % nm
rmax=rPFunc(rmaxnm*10^-9);        % dimensionless "rbar"
nr=         700;                % number of r spatial positions 
nt=         1000;               % number of time incrments
freqHz = 1e0;                    % Hz ... frequency of electron pulses
freq = freqHz/tPFunc(1);        % dimensionless ... frequency of electron pulses
duty = 1;                     % duty cycle of electron pulses

%---window
range = rmax;                  % position plot range
edge = -rmax;                  % position plot offset

%---pulse parameters

fprintf('\n');
fprintf('Simulation Parameters\n');
fprintf('T                     = %g\n',T);
fprintf('nr                    = %g\n',nr);
fprintf('nt                    = %g\n',nt);
fprintf('rmax                  = %g\n',rmax);
fprintf('\n');