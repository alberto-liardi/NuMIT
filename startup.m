%% Startup script for the NUMIT package
%
% Initialise package and add subfolder dependencies to the path.
% This script is run automatically if MATLAB is launched in the package directory.

disp("Initialising NuMIT package dependencies...")
parts = fileparts(mfilename('fullpath')); % directory containing this file

addpath(parts+"/core");
addpath(parts+"/io");
addpath(parts+"/utils");
addpath(parts+"/private/MVGC2");
addpath(parts+"/private/partial-info-decomp-master");

% now initialising MVGC2
startupMVGC2 = parts+"/private/MVGC2/startup.m";
run(startupMVGC2);

disp("Initialisation complete!");