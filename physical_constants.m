%---physical constants
global T ii g G D B_r c T1 rmax u0 T12 T13 mu_p mu_e B0 B_d kB temp freq duty

ge = -2.00231930436153;
gp = 5.585694713;
hb = 1.054571726e-34; % m^2 kg/sec
gamma_e = 1.760859708e11; % nuclear gyromagnetic ratio, T^(-1) s^(-1)
gamma_p = 2.675222005e8; % electron gyromagnetic ratio, T^(-1) s^(-1)
mu = 4*pi*1e-7; % T m/A ... magnetic constant
kB = 1.3806488e-23; % J/K ... Boltzmann constant
NA = 6.02214129e23; % mol^-1 ... Avogadro constant

mu_B = -hb*gamma_e/ge; % Bohr magneton
mu_e = -hb*gamma_e/2; % mag moment of electron spin
mu_N = hb*gamma_p/gp; % nuclear magneton
mu_p = hb*gamma_p/2; % nuclear spin magnetic moment
