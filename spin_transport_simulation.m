classdef spin_transport_simulation
   properties
      class_version = 0.1; % float for easy version comparison
      class_definition % to save a copy of this classdef file for comparison of serialized instances
      timestamp = datestr(now); % timestamp of class creation
      physical_constants
      pde
      initial_conditions = {};
      boundary_conditions = {};
      ode_solver_options = {};
      simulation_parameters = {};
      simulation_results = {};
   end
   methods
       function obj = spin_transport_simulation() % called at instance creation
             obj.class_definition = fileread([mfilename(),'.m']);
       end
   end
end