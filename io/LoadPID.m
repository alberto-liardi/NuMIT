%% LoadPID
% 
% Loads PID atoms from a .csv file and computes the total mutual information (MI).
% 
% Description:
%   This function reads PID atoms from a specified file and returns arrays for each of the 
%   PID atoms: Unique X, Unique Y, Redundancy, and Synergy. 
%   It also computes the total mutual information (MI) for each sample.
% 
% Inputs:
%   file    - Path to the file containing PID atom data in table format.
% 
% Outputs:
%   UnX     - [1xL] array where each element is the Unique to X atom for the corresponding sample.
%   UnY     - [1xL] array where each element is the Unique to Y atom for the corresponding sample.
%   Red     - [1xL] array where each element is the Redundancy atom for the corresponding sample.
%   Syn     - [1xL] array where each element is the Synergy atom for the corresponding sample.
%   MI      - [1xL] array where each element is the total mutual information for the corresponding sample.
% 
% Alberto Liardi, 2024


function [UnX, UnY, Red, Syn, MI] = LoadPID(file)

    if ~isfile(file), error(file+" directory has not been found!"); end

    T = readtable(file);
    PIDs = table2array(T);
    UnX = PIDs(:,1)'; UnY = PIDs(:,2)'; Red = PIDs(:,3)'; Syn = PIDs(:,4)';
    MI = sum(PIDs,2)';

end