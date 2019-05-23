classdef spin_transport_simulation < handle % enables self-updating
  properties
    class_version = 0.1; % float for easy version comparison
    class_definition % to save a copy of this classdef file for comparison of serialized instances
    timestamp = datestr(now); % timestamp of class creation
    physical_constants
    pde
    initial_conditions = @(self,rr) initial_conditions_equilibrium(self,rr); % set initial conditions method
    boundary_conditions = {};
    ode_solver_options = {};
    simulation_parameters = {};
    simulation_results = {};
  end
  methods
    function self = spin_transport_simulation() % called at instance creation
      self.class_definition = fileread([mfilename(),'.m']);
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
  end
end