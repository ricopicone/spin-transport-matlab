classdef spin_transport_simulation < handle % enables self-updating
  properties
    class_version = 0.1; % float for easy version comparison
    class_definition % to save a copy of this classdef file for comparison of serialized instances
    timestamp = datestr(now); % timestamp of class creation
    constants = struct(); % physical constants
    constants_set = @(self) constants_nominal(self); % set physical constants
    pde = @(self,r,t,p,DpDx) pde_nominal(self,r,t,p,DpDx) % pde to use in solution
    initial_conditions = @(self,rr) initial_conditions_equilibrium(self,rr); % set initial conditions method
    boundary_conditions = @(self,xl,ul,xr,ur,t) boundary_conditions_equilibrium_j(self,xl,ul,xr,ur,t); % set boundary conditions method
    ode_solver_options = struct('AbsTol',1e-12,'RelTol',1e-6);
    parameters = struct();
    parameters_set = @(self,rr) parameters_nominal(self); % set initial conditions method
    results = struct();
    docs = struct(); % Documentation of "sub-properties" of the class.
  end
  methods
    function self = spin_transport_simulation() % called at instance creation
      self.docs.constants = struct(); % initialize
      self.docs.parameters = struct(); % initialize
      self.class_definition = fileread([mfilename(),'.m']);
      self.constants_set(self); % define constants
      self.parameters_set(self); % define parameters
    end
    function self = constants_nominal(self)
      c = self.constants; % unpack
      d = self.docs.constants; % unpack
      c.ge = -2.00231930436153;
      d.ge = [...
        'electron spin g-factor\n'...
        'units: dimensionless\n'...
        'source: https://en.wikipedia.org/wiki/G-factor_(physics)#Electron_spin_g-factor'...
      ];
      c.gp = 5.585694713;
      d.gp = [...
        'nuclear spin g-factor\n'...
        'units: dimensionless\n'...
        'source: http://leona.physics.tamu.edu/Phys327.11s/lab-gyro.pdf'...
      ];
      c.hb = 1.054571726e-34;
      d.hb = [...
        'reduced Plancks constant\n'...
        'units: m^2 kg/s\n'...
        'source: https://en.wikipedia.org/wiki/Planck_constant#Value'
      ];
      c.gamma_e = 1.760859708e11; % nuclear gyromagnetic ratio, T^(-1) s^(-1)
      c.gamma_p = 2.675222005e8; % electron gyromagnetic ratio, T^(-1) s^(-1)
      c.mu = 4*pi*1e-7; % T m/A ... magnetic constant
      c.kB = 1.3806488e-23; % J/K ... Boltzmann constant
      c.NA = 6.02214129e23; % mol^-1 ... Avogadro constant
      c.mu_B = -c.hb*c.gamma_e/c.ge; % Bohr magneton
      c.mu_e = -c.hb*c.gamma_e/2; % mag moment of electron spin
      c.mu_N = c.hb*c.gamma_p/c.gp; % nuclear magneton
      c.mu_p = c.hb*c.gamma_p/2; % nuclear spin magnetic moment
      self.constants = c; % repack
      self.docs.constants = d; % repack
    end
    function self = parameters_nominal(self)
      p = self.parameters; % unpack
      c = self.constants; % unpack
      d = self.docs.parameters; % unpack
      % system parameters
      p.MwPS = 12*0 + 1*8;  %  no spin-1/2 from C atoms; H atoms
                          %  g/mol ... molar mass of polystyrene assuming C8H8 (Wikipedia)
      p.dPS=1.047; % g/mL ... density of polystyrene (from bottle)
      p.nAMPS = 2; % number of ^1H atoms per molecule for polystyrene
      p.concDPPH = .01; % concentration DPPH
      p.concPS = 1-p.concDPPH; % concentration polystyrene *)
      p.den2 = 1/p.MwPS*p.dPS*1e6*c.NA*p.nAMPS; % 1/m^3
      p.Delta_2 = p.concPS*p.den2; % 1/m^3  \[Delta]2 -> 1.34*10^26 1/m^3, Dougherty2000 
      p.MwDPPH = 394.32; % g/mol ... molar mass of DPPH (wikipedia)
      p.dDPPH = 1.4; % g/cm^3 ... density of DPPH (wikipedia)
      p.nAMDPPH = 1; % just one free radical per molecule
      p.den3 = 1/p.MwDPPH*p.dDPPH*1e6*c.NA*p.nAMDPPH; % 1/m^3
      p.Delta_3 = p.concDPPH*p.den3; % 1/m^3
      p.Gamma_2 = c.mu/(4*pi)*c.hb*c.gamma_p^2*p.Delta_2^(1/3); % rad/(sec Tesla)
      p.Gamma_3 = c.mu/(4*pi)*c.hb*c.gamma_e^2*p.Delta_3^(1/3); % rad/(sec Tesla)
      p.grad = 1e2*44000; % magnetic field gradient (T/m);
      p.Bd_2 = c.mu/(4*pi)*c.hb*c.gamma_p*p.Delta_2;    % T 
      p.Bd_3 = c.mu/(4*pi)*c.hb*c.gamma_e*p.Delta_3;    % T 
      p.B_d = (p.Bd_2 + p.Bd_3);  % Bd2 + Bd3
      p.B0 = 2.7; % T
      p.B1max_p_nom = 1e-3; % T ... TODO I made this up
      p.B1max_e_nom = 1e-3; % T ... TODO I made this up
      p.temp = 10; % K
      p.tPFunc = @(t) p.Gamma_2*(p.grad/p.B_d)^2*t;
      p.rPFunc = @(r) p.grad/p.B_d*r;
      p.T12sec = 0.1; % sec ... T1 proton
      p.T13sec = 30.3 * 10^-6; % sec ... T1 electron
      p.T12 = p.tPFunc(p.T12sec);
      p.T13 = p.tPFunc(p.T13sec);
      p.T22sec = 20e-9; % sec ... T2 proton ... TODO made this up
      p.T23sec = 20e-9; % sec ... T2 electron ... TODO made this up
      p.T22 = p.tPFunc(p.T22sec);
      p.T23 = p.tPFunc(p.T23sec);
      % dimensionless system parameters
      p.g = c.gamma_e/c.gamma_p;    
      p.G = p.Gamma_3/p.Gamma_2;    
      p.D = p.Delta_3/p.Delta_2;    
      p.B_r = 1;    
      p.c = p.B_r*(1+p.D)/(1+p.g*p.D);
      % simulation parameters
      p.t_max_sec = 2e-9*p.G;
      p.t_max = p.tPFunc(p.t_max_sec); % dimensionless ... normalized max simulation time
      p.r_max_nm = 2; % nm
      p.r_max = p.rPFunc(p.r_max_nm*1e-9); % dimensionless "rbar"
      p.n_r = 400; % number of r spatial positions
      p.n_r = 1000; % number of time increments
      p.pulse.n_p = 1; % number of proton pulses
      p.pulse.n_e = 1; % number of electron pulses
      p.pulse.duty_p = .1; % duty cycle of proton pulses
      p.pulse.duty_e = .01; % duty cycle of electron pulses
      p.pulse.on_B1_p = 0; % nuclear B1 on?
      p.pulse.on_B1_e = 0; % electron B1 on?
      p.pulse.B1_max_p = p.pulse.on_B1_p*p.B1max_p_nom;
      p.pulse.B1_max_e = p.pulse.on_B1_e*p.B1max_e_nom;
      p.pulse.T_p = p.t_max/p.pulse.n_p; % dimensionless ... period of proton pulses
      p.pulse.freq_p = 1/p.pulse.T_p; % dimensionless ... frequency of proton pulses
      p.pulse.T_s_p_nice = 1/(4*p.pulse.freq_p); % nice period of sinusoidal transition for B1 pulse
      p.pulse.T_s_p = min([... % period of sinusoidal transition for B1 pulse
        p.pulse.T_s_p_nice,...
        2*p.pulse.duty_p/p.pulse.freq_p,...
        2*(1-p.pulse.duty_p)/p.pulse.freq_p...
      ]); 
      p.pulse.T_e = p.t_max/p.pulse.n_e; % dimensionless ... period of electron pulses
      p.pulse.freq_e = 1/p.pulse.T_e; % dimensionless ... frequency of electron pulses
      p.pulse.T_s_e_nice = 1/(4*p.pulse.freq_e); % nice period of sinusoidal transition for B1 pulse
      p.pulse.T_s_e = min([... % period of sinusoidal transition for B1 pulse
        p.pulse.T_s_e_nice,...
        2*p.pulse.duty_e/p.pulse.freq_e,...
        2*(1-p.pulse.duty_e)/p.pulse.freq_e...
      ]); 
      p.plot.range = 1*p.r_max/2; % position plot range
      p.plot.edge = -1*p.r_max/2; % position plot offset
      % langevin
      mu_p_kB = c.mu_p/c.kB;
      mu_e_kB = c.mu_e/c.kB;
      p.rho_1_langevin= @(r) ... % langevin energy
        -1*... % I think this is because mag grad is negative
          tanh(...
            (...
                1 + c.gamma_e/c.gamma_p * p.Delta_3/p.Delta_2 ...
            )/(c.gamma_e/c.gamma_p * ( 1 + p.Delta_3/p.Delta_2 ) ) * ...
            (...
              c.mu_e * p.B_d / ...
              c.kB * p.temp ...
            ) ...
          );
      p.rho_2_langevin= @(r) tanh((p.B0-p.B_d.*r).*mu_p_kB./p.temp); % langevin nuclear polarization
      p.rho_3_langevin= @(r) tanh((p.B0-p.B_d.*r).*mu_e_kB./p.temp); % langevin electron polarization
      % gradient of equilibrium solution (langevin)
      p.d_rho_1_langevin= @(r) 0;
      p.d_rho_2_langevin= @(r) double(-p.B_d*mu_p_kB/p.temp*(sech((p.B0-p.B_d*r)*mu_p_kB/p.temp))^2);
      p.d_rho_3_langevin= @(r) double(-p.B_d*mu_e_kB/p.temp*(sech((p.B0-p.B_d*r)*mu_e_kB/p.temp))^2);
      self.parameters = p; % repack
      self.docs.parameters = d; % repack
    end
    function u0 = initial_conditions_equilibrium(self,rr)
      % just zeros!
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
      G = self.parameters.G;
      ql = [-1/(1+G);-1;-1/G];                                                         
      qr = [-1/(1+G);-1;-1/G];  
      pl = [d_rho_1_langevin(xl);d_rho_2_langevin(xl);d_rho_3_langevin(xl)];
      pr = [d_rho_1_langevin(xr);d_rho_2_langevin(xr);d_rho_3_langevin(xr)];
    end
    function [pl,ql,pr,qr] = boundary_conditions_zero_j(self,xl,ul,xr,ur,t)
      G = self.parameters.G;
      c = self.parameters.c;
      g = self.parameters.g;
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
    function [cm,f,s] = pde_nominal(self,r,t,p,DpDx)
      c = self.constants; % unpack
      p = self.parameters; % unpack
      cm = [1; 1; 1];                                  
      f = [...
        1+p.G; ...
        1; ...
        p.G ...
      ].*DpDx;

      % separation terms
      % s = [...
      %     -c^2/(1+D)*( (1-p(2)^2) + G*D*g^2*(1-p(3)^2) )*atanh(p(1))-c/(1+D)*(DpDx(2)-G*D*g*DpDx(3)); ...
      %     c*( (1-p(2)^2)/(1-p(1)^2)*DpDx(1) - 2*atanh(p(1))*p(2)*DpDx(2) ); ...
      %     -c*G*g*( (1-p(3)^2)/(1-p(1)^2)*DpDx(1) - 2*atanh(p(1))*p(3)*DpDx(3) ) ...
      % ];  
      s = [0;0;0]; % no cross terms

      % Bloch dynamics sans relaxation
      % - amplitude modulation
      % w1_2 = -gamma_p*B1(t,freq_p,duty_p,T_s_p,B1max_p);
      % w1_3 = -gamma_e*B1(t,freq_e,duty_e,T_s_e,B1max_e);
      % - frequency modulation
      % % delta_tilde_2 = gamma_p*Bd_2*(r-rs);
      % % delta_tilde_3 = gamma_e*Bd_3*(r-rs);
      % delta_tilde_2 = gamma_p*Bd_2*r;
      % delta_tilde_3 = gamma_e*Bd_3*r;
      % - modulation enters through the taus
      % tau_2 = 1/(1/T12+T22*w1_2^2/(1+T22^2*delta_tilde_2^2));
      % tau_3 = 1/(1/T13+T23*w1_3^2/(1+T23^2*delta_tilde_3^2));
      % s = s + [ ...
      %     0; ...
      %     -p(2)/tau_2; ...
      %     -p(3)/tau_3 ...
      % ];

      % Bloch relaxation dynamics <--this slows the solver
      % s = s + [ ...
      %     0; ...
      %     rho_2_langevin(r)/T12; ...
      %     rho_3_langevin(r)/T13 ...
      % ];

      % print
      ii = ii+1;
      if mod(ii,10000)==0
          fprintf('%d  %g\n',ii,t)
      end
    end
    function [cm,f,s] = pde_diffusion(self,r,t,p,DpDx)
    % PDE_DIFFUSION  pde for diffusion unit test
      c = self.constants; % unpack
      p = self.parameters; % unpack
      cm = [1; 1; 1];                                  
      f = [...
        1+p.G; ...
        1; ...
        p.G ...
      ].*DpDx;
      s = [0;0;0]; % no cross terms
    end
  end
end