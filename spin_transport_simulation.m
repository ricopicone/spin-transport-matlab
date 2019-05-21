classdef spin_transport_simulation
   properties
      timestamp
      physical_constants
      pde
      initial_conditions
      boundary_conditions
      ode_solver_options
      simulation_parameters
   end
   methods
       function obj = spin_transport_simulation()
          if nargin == 0
             obj.timestamp = datestr(now);
          end
       end
   end
end