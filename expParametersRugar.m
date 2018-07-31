global rho_1_langevin rho_2_langevin rho_3_langevin
%---experimental parameters
% Delta_2=    1.51324e29;
% Delta_3=    8.55244e25;
% Gamma_2=    4.02191e-15;        % rad/(sec Tesla)
% Gamma_3=    1.44064e-10;        % rad/(sec Tesla)
MwPS = 12*0 + 1*8;  %  no spin-1/2 from C atoms; H atoms
                    %  g/mol ... molar mass of polystyrene assuming C8H8 (Wikipedia)
dPS=1.047;          % g/mL ... density of polystyrene (from bottle)
nAMPS = 2;          % number of ^1H atoms per molecule for polystyrene
concDPPH = .01;     % concentration DPPH
concPS = 1-concDPPH; % concentration polystyrene *)
den2 = 1/MwPS*dPS*1e6*NA*nAMPS; % 1/m^3
Delta_2 = concPS*den2; % 1/m^3  \[Delta]2 -> 1.34*10^26 1/m^3, Dougherty2000 
MwDPPH = 394.32;    % g/mol ... molar mass of DPPH (wikipedia)
dDPPH = 1.4;        % g/cm^3 ... density of DPPH (wikipedia)
nAMDPPH = 1;        % just one free radical per molecule
den3 = 1/MwDPPH*dDPPH*1e6*NA*nAMDPPH; % 1/m^3
Delta_3 = concDPPH*den3; % 1/m^3
Gamma_2=mu/(4*pi)*hb*gamma_p^2*Delta_2^(1/3);        % rad/(sec Tesla)
Gamma_3=mu/(4*pi)*hb*gamma_e^2*Delta_3^(1/3);        % rad/(sec Tesla)
grad=       1e2*44000;              % magnetic field gradient (T/m);
Bd_2 = mu/(4*pi)*hb*gamma_p*Delta_2;    % T 
Bd_3 = mu/(4*pi)*hb*gamma_e*Delta_3;    % T 
B_d = (Bd_2 + Bd_3);  % Bd2 + Bd3
%    B_d=        0.0000585733;       % T
% B0=         0.0893;             % t
B0=         2.7;             % t
temp=       10;                 % K
tPFunc = @(t) Gamma_2*(grad/B_d)^2*t;
rPFunc = @(r) grad/B_d * r;
T12sec = 0.1;           % sec ... T1 proton
T13sec = 30.3 * 10^-6;  % sec ... T1 electron
T12 = tPFunc(T12sec);
T13 = tPFunc(T13sec);

%---Langevin
rho_1_langevin= @(r)  ( 1 + gamma_e/gamma_p * Delta_3/Delta_2 ) * mu_p * B_d / ( kB * temp );    % langevin nuclear polarization
rho_2_langevin= @(r)  tanh((B0-B_d*r)*mu_p/kB/temp);    % langevin nuclear polarization
rho_3_langevin= @(r)  tanh((B0-B_d*r)*mu_e/kB/temp);    % langevin electron polarization

fprintf('\n');
fprintf('Experimental Parameters\n');
fprintf('concDPPH              = %g\n',concDPPH);
fprintf('Delta_2               = %g\n',Delta_2);
fprintf('Delta_3               = %g\n',Delta_3);
fprintf('Gamma_2               = %grad/(sec Tesla)\n',Gamma_2);
fprintf('Gamma_3               = %g rad/(sec Tesla)\n',Gamma_3);
fprintf('B_d                   = %g T\n',B_d);
fprintf('B0                    = %g T\n',B0);
fprintf('temp                  = %g K\n',temp);