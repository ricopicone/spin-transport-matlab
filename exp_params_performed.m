global rho_1_langevin rho_2_langevin rho_3_langevin
%---experimental parameters
% Delta_2=    1.51324e29;
% Delta_3=    8.55244e25;
% Gamma_2=    4.02191e-15;        % rad/(sec Tesla)
% Gamma_3=    1.44064e-10;        % rad/(sec Tesla)

% general

temp = 10; % K ... physical temp of sample

% electron spin concentration
MwDPPH = 394.32;    % g/mol ... molar mass of DPPH (wikipedia)
dDPPH = 1.4;        % g/cm^3 ... density of DPPH (wikipedia)
nAMDPPH = 1;        % just one free radical per molecule
concDPPH = .04;     % concentration DPPH
den3 = 1/MwDPPH*dDPPH*1e6*NA*nAMDPPH; % 1/m^3
Delta_3 = concDPPH*den3; % 1/m^3

% nuclear spin concentration
MwPS = ... % g/mol ... molar mass of polystyrene assuming C8H8 (Wikipedia)
  12*0 + ... % no spin-1/2 from C atoms
  1*8; % H atoms
dPS=1.047; % g/mL ... density of polystyrene (from bottle)
nAMPS = 2; % number of ^1H atoms per molecule for polystyrene
concPS = 1-concDPPH; % concentration polystyrene
den2 = 1/MwPS*dPS*1e6*NA*nAMPS; % 1/m^3
Delta_2 = concPS*den2; % 1/m^3 ... Delta_2 = 1.34*10^26 1/m^3, Dougherty2000 

% identifying spins
gamma_2 = gamma_p; % species 2 is proton
gamma_3 = gamma_e; % species 3 is electron

% magnetic fields
B0 = 0.0893; % T ... background field
grad = 44e3; % T/m ... magnetic field gradient;
Bd_2 = mu/(4*pi)*hb*gamma_2*Delta_2; % T ... nuclear dipole field
Bd_3 = mu/(4*pi)*hb*gamma_3*Delta_3; % T ... electron dipole field
B_d = Bd_2 + Bd_3;  % T ... not the worst guess

% transport rates
Gamma_2=mu/(4*pi)*hb*gamma_2^2*Delta_2^(1/3); % rad/(sec Tesla)
Gamma_3=mu/(4*pi)*hb*gamma_3^2*Delta_3^(1/3); % rad/(sec Tesla)

% relaxation rates
T1e = 30.3e-6;  % sec ... electron T1 at 10K (Dougherty2000)
T2e = 20e-9;    % sec ... electron T2 at 10K (Dougherty2000)
T1p = .1;       % sec ... nuclear T1
T2p = T2e;      % sec ... nuclear T2 ... TODO look up good value
T12 = T1p
T22 = T2p
T13 = T1e
T23 = T2e


% print
% fprintf('\n');
% fprintf('Experimental Parameters\n');
% fprintf('concDPPH              = %g\n',concDPPH);
% fprintf('Delta_2               = %g\n',Delta_2);
% fprintf('Delta_3               = %g\n',Delta_3);
% fprintf('Gamma_2               = %grad/(sec Tesla)\n',Gamma_2);
% fprintf('Gamma_3               = %g rad/(sec Tesla)\n',Gamma_3);
% fprintf('B_d                   = %g T\n',B_d);
% fprintf('B0                    = %g T\n',B0);
% fprintf('temp                  = %g K\n',temp);