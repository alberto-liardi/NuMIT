%% PID_DEP_Gaussian
% 
% Compute Partial Information Decomposition (PID) components for a Gaussian system using
% the Unique Information via Dependency Constraines (DEP) definition.
% 
% Description:
%   This function calculates the DEP-PID atoms between two sources and a target 
%   in a Gaussian system. Sources and target variables can be multivariate.
% 
% Inputs:
%   Sigma_full - [S+T, S+T] covariance matrix of the full system.
%   S          - Number of sources in the system.
%   T          - Number of targets in the system.
% 
% Outputs:
%   UnX        - Unique information from the first source.
%   UnY        - Unique information from the second source.
%   Red        - Redundancy.
%   Syn        - Synergy.
% 
% Alberto Liardi, 2024

function [UnX, UnY, Red, Syn] = PID_DEP_Gaussian(Sigma_full,S,T)

    PI = num2cell(calc_pi_Idep_mvn(Sigma_full, [S/2 S/2 T])/log2(exp(1)));
    [Red, UnX, UnY, Syn] = PI{:};

end