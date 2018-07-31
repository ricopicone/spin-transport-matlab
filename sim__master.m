clear
clear global
close all;

global sim

%% Choose simulation switch

% sim = 'performed';  % experiment performed
sim = 'rugar1';     % using Rugar's field B0 and gradient

save_flag = 0;

%% Physical constants

physical_constants;

%% Experimental parameters

switch sim
    case 'performed'
        exp_params_performed;
    case 'rugar1'
        exp_params_rugar;
end

%% Normalized parameters

normalized_params;

%% Simulation parameters

switch sim
    case 'performed'
        sim_params_performed;
    case 'rugar1'
        sim_params_rugar;
end

%% Grid

grid_compute;

%% Simulate

n_traces = 20;
sim_simulate;

%% Postprocess

postprocess;

%% Plots and results

plots_results;