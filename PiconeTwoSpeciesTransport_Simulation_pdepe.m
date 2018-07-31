clear
clear global
close all;

global sim

%% Choose simulation switch
% sim = 'performed';  % experiment performed
sim = 'rugar1';     % using Rugar's field B0 and gradient

saveFlag = 0;

%% Physical constants
physicalConstants;

%% Experimental parameters
switch sim
    case 'performed'
        expParametersPerformedExperiment;
    case 'rugar1'
        expParametersRugar;
end

%% Normalized parameters
normalizedParameters;

%% Simulation parameters
switch sim
    case 'performed'
        simParametersPerformedExperiment;
    case 'rugar1'
        simParametersRugar1;
end

%% Grid
gridCompute;

%% Simulate
nTraces = 20;
simulate;

%% Plots and results
plotsAndResults;