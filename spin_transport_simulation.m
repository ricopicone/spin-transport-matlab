classdef spin_transport_simulation < handle % enables self-updating
  properties
    class_version = 0.1; % float for easy version comparison
    class_definition % to save a copy of this classdef file for comparison of serialized instances
    timestamp = datestr(now); % timestamp of class creation
    constants = struct(); % physical constants
    constants_set = @(self) constants_nominal(self); % set physical constants
    pde
    initial_conditions = @(self,rr) initial_conditions_equilibrium(self,rr); % set initial conditions method
    boundary_conditions = @(self,xl,ul,xr,ur,t) boundary_conditions_equilibrium_j(self,xl,ul,xr,ur,t); % set boundary conditions method
    ode_solver_options = struct( ...
      'AbsTol',1e-12, ...
      'RelTol',1e-6 ...
    );
    parameters = struct();
    parameters_set = @(self,rr) parameters_nominal(self); % set initial conditions method
    results = struct();
  end
  methods
    function self = spin_transport_simulation() % called at instance creation
      self.class_definition = fileread([mfilename(),'.m']);
      self.constants_set(self); 
      self.parameters_set(self);
    end
    function self = constants_nominal(self)
      self.constants.ge = -2.00231930436153;
      self.constants.gp = 5.585694713;
      self.constants.hb = 1.054571726e-34; % m^2 kg/sec
      self.constants.gamma_e = 1.760859708e11; % nuclear gyromagnetic ratio, T^(-1) s^(-1)
      self.constants.gamma_p = 2.675222005e8; % electron gyromagnetic ratio, T^(-1) s^(-1)
      self.constants.mu = 4*pi*1e-7; % T m/A ... magnetic constant
      self.constants.kB = 1.3806488e-23; % J/K ... Boltzmann constant
      self.constants.NA = 6.02214129e23; % mol^-1 ... Avogadro constant
      self.constants.mu_B = -self.constants.hb*self.constants.gamma_e/self.constants.ge; % Bohr magneton
      self.constants.mu_e = -self.constants.hb*self.constants.gamma_e/2; % mag moment of electron spin
      self.constants.mu_N = self.constants.hb*self.constants.gamma_p/self.constants.gp; % nuclear magneton
      self.constants.mu_p = self.constants.hb*self.constants.gamma_p/2; % nuclear spin magnetic moment
    end
    function self = parameters_nominal(self)
      % system parameters
      self.parameters.MwPS = 12*0 + 1*8;  %  no spin-1/2 from C atoms; H atoms
                          %  g/mol ... molar mass of polystyrene assuming C8H8 (Wikipedia)
      self.parameters.dPS=1.047; % g/mL ... density of polystyrene (from bottle)
      self.parameters.nAMPS = 2; % number of ^1H atoms per molecule for polystyrene
      self.parameters.concDPPH = .01; % concentration DPPH
      self.parameters.concPS = 1-self.parameters.concDPPH; % concentration polystyrene *)
      self.parameters.den2 = 1/self.parameters.MwPS*self.parameters.dPS*1e6*self.constants.NA*self.parameters.nAMPS; % 1/m^3
      self.parameters.Delta_2 = self.parameters.concPS*self.parameters.den2; % 1/m^3  \[Delta]2 -> 1.34*10^26 1/m^3, Dougherty2000 
      self.parameters.MwDPPH = 394.32; % g/mol ... molar mass of DPPH (wikipedia)
      self.parameters.dDPPH = 1.4; % g/cm^3 ... density of DPPH (wikipedia)
      self.parameters.nAMDPPH = 1; % just one free radical per molecule
      self.parameters.den3 = 1/self.parameters.MwDPPH*self.parameters.dDPPH*1e6*self.constants.NA*self.parameters.nAMDPPH; % 1/m^3
      self.parameters.Delta_3 = self.parameters.concDPPH*self.parameters.den3; % 1/m^3
      self.parameters.Gamma_2 = self.constants.mu/(4*pi)*self.constants.hb*self.constants.gamma_p^2*self.parameters.Delta_2^(1/3); % rad/(sec Tesla)
      self.parameters.Gamma_3 = self.constants.mu/(4*pi)*self.constants.hb*self.constants.gamma_e^2*self.parameters.Delta_3^(1/3); % rad/(sec Tesla)
      self.parameters.grad = 1e2*44000; % magnetic field gradient (T/m);
      self.parameters.Bd_2 = self.constants.mu/(4*pi)*self.constants.hb*self.constants.gamma_p*self.parameters.Delta_2;    % T 
      self.parameters.Bd_3 = self.constants.mu/(4*pi)*self.constants.hb*self.constants.gamma_e*self.parameters.Delta_3;    % T 
      self.parameters.B_d = (self.parameters.Bd_2 + self.parameters.Bd_3);  % Bd2 + Bd3
      self.parameters.B0 = 2.7; % T
      self.parameters.B1max_p_nom = 1e-3; % T ... TODO I made this up
      self.parameters.B1max_e_nom = 1e-3; % T ... TODO I made this up
      self.parameters.temp = 10; % K
      self.parameters.tPFunc = @(t) self.parameters.Gamma_2*(self.parameters.grad/self.parameters.B_d)^2*t;
      self.parameters.rPFunc = @(r) self.parameters.grad/self.parameters.B_d*r;
      self.parameters.T12sec = 0.1; % sec ... T1 proton
      self.parameters.T13sec = 30.3 * 10^-6; % sec ... T1 electron
      self.parameters.T12 = self.parameters.tPFunc(self.parameters.T12sec);
      self.parameters.T13 = self.parameters.tPFunc(self.parameters.T13sec);
      self.parameters.T22sec = 20e-9; % sec ... T2 proton ... TODO made this up
      self.parameters.T23sec = 20e-9; % sec ... T2 electron ... TODO made this up
      self.parameters.T22 = self.parameters.tPFunc(self.parameters.T22sec);
      self.parameters.T23 = self.parameters.tPFunc(self.parameters.T23sec);
      % dimensionless system parameters
      self.parameters.g = self.constants.gamma_e/self.constants.gamma_p;    
      self.parameters.G = self.parameters.Gamma_3/self.parameters.Gamma_2;    
      self.parameters.D = self.parameters.Delta_3/self.parameters.Delta_2;    
      self.parameters.B_r = 1;    
      self.parameters.c = self.parameters.B_r*(1+self.parameters.D)/(1+self.parameters.g*self.parameters.D);
      % simulation parameters
      self.parameters.t_max_sec = 2e-9*self.parameters.G;
      self.parameters.t_max = self.parameters.tPFunc(self.parameters.t_max_sec); % dimensionless ... normalized max simulation time
      self.parameters.r_max_nm = 2; % nm
      self.parameters.r_max = self.parameters.rPFunc(self.parameters.r_max_nm*1e-9); % dimensionless "rbar"
      self.parameters.n_r = 400; % number of r spatial positions
      self.parameters.n_r = 1000; % number of time incrments
      self.parameters.pulse.n_p = 1; % number of proton pulses
      self.parameters.pulse.n_e = 1; % number of electron pulses
      self.parameters.pulse.duty_p = .1; % duty cycle of proton pulses
      self.parameters.pulse.duty_e = .01; % duty cycle of electron pulses
      self.parameters.pulse.on_B1_p = 0; % nuclear B1 on?
      self.parameters.pulse.on_B1_e = 0; % electron B1 on?
      self.parameters.pulse.B1_max_p = self.parameters.pulse.on_B1_p*self.parameters.B1max_p_nom;
      self.parameters.pulse.B1_max_e = self.parameters.pulse.on_B1_e*self.parameters.B1max_e_nom;
      self.parameters.pulse.T_p = self.parameters.t_max/self.parameters.pulse.n_p; % dimensionless ... period of proton pulses
      self.parameters.pulse.freq_p = 1/self.parameters.pulse.T_p; % dimensionless ... frequency of proton pulses
      self.parameters.pulse.T_s_p_nice = 1/(4*self.parameters.pulse.freq_p); % nice period of sinusoidal transition for B1 pulse
      self.parameters.pulse.T_s_p = min([... % period of sinusoidal transition for B1 pulse
        self.parameters.pulse.T_s_p_nice,...
        2*self.parameters.pulse.duty_p/self.parameters.pulse.freq_p,...
        2*(1-self.parameters.pulse.duty_p)/self.parameters.pulse.freq_p...
      ]); 
      self.parameters.pulse.T_e = self.parameters.t_max/self.parameters.pulse.n_e; % dimensionless ... period of electron pulses
      self.parameters.pulse.freq_e = 1/self.parameters.pulse.T_e; % dimensionless ... frequency of electron pulses
      self.parameters.pulse.T_s_e_nice = 1/(4*self.parameters.pulse.freq_e); % nice period of sinusoidal transition for B1 pulse
      self.parameters.pulse.T_s_e = min([... % period of sinusoidal transition for B1 pulse
        self.parameters.pulse.T_s_e_nice,...
        2*self.parameters.pulse.duty_e/self.parameters.pulse.freq_e,...
        2*(1-self.parameters.pulse.duty_e)/self.parameters.pulse.freq_e...
      ]); 
      self.parameters.plot.range = 1*self.parameters.r_max/2; % position plot range
      self.parameters.plot.edge = -1*self.parameters.r_max/2; % position plot offset
    end
    function u0 = initial_conditions_equilibrium(self,rr)
      u0 = [0;0;0];
    end
    function u0 = initial_conditions_guassian(self,rr)
      u0 = ...
        [...
          0; ...
          1 * ( exp( - ( rr )^2/1e-1 ) ); ...
          1 * ( exp( - ( rr )^2/1e-1 ) ) ...
        ];
    end
    function [pl,ql,pr,qr] = boundary_conditions_equilibrium_j(self,xl,ul,xr,ur,t)
      ql = [-1/(1+G);-1;-1/G];                                                         
      qr = [-1/(1+G);-1;-1/G];  
      pl = [d_rho_1_langevin(xl);d_rho_2_langevin(xl);d_rho_3_langevin(xl)];
      pr = [d_rho_1_langevin(xr);d_rho_2_langevin(xr);d_rho_3_langevin(xr)];
    end
    function [pl,ql,pr,qr] = boundary_conditions_zero_j(self,xl,ul,xr,ur,t)
      pl = [0; -c*(1-ul(2)^2)*atanh(ul(1)); -c*G*g*(1-ul(3)^2)*atanh(ul(1))];                               
      ql = [-1;-1;-1];                                  
      pr = [0; -c*(1-ur(2)^2)*atanh(ur(1)); -c*G*g*(1-ur(3)^2)*atanh(ur(1))];                         
      qr = [-1;-1;-1]; 
    end
    function [pl,ql,pr,qr] = boundary_conditions_zero_gradient(self,xl,ul,xr,ur,t)
      pl = [0; 0; 0];                               
      ql = [1;1;1];                                  
      pr = [0; 0; 0];                         
      qr = [1;1;1]; 
    end
    function [pl,ql,pr,qr] = boundary_conditions_pinned(self,xl,ul,xr,ur,t)
      pl = ul - [rho_1_langevin(xl);rho_2_langevin(xl);rho_3_langevin(xl)];
      ql = [0;0;0];
      pr = ur - [rho_1_langevin(xr);rho_2_langevin(xr);rho_3_langevin(xr)];
      qr = [0;0;0];  
    end
  end
end