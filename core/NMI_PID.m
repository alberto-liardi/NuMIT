%% NMI_PID
% 
% Normalises PID atoms with respect to the corresponding mutual information.
% 
% Description:
%   This function normalizes each of the PID atoms (Unique X, Unique Y, Redundancy, Synergy)
%   by dividing them by the corresponding mutual information (MI) value. 
% 
% Inputs:
%   Un_X  - Vector of Unique X PID atoms.
%   Un_Y  - Vector of Unique Y PID atoms.
%   Red   - Vector of Redundancy PID atoms.
%   Syn   - Vector of Synergy PID atoms.
%   MI    - Vector of mutual information values.
% 
% Outputs:
%   U_X   - Vector of normalised Unique X PID atoms.
%   U_Y   - Vector of normalised Unique Y PID atoms.
%   R     - Vector of normalised Redundancy PID atoms.
%   S     - Vector of normalised Synergy PID atoms.
%   MI    - Vector of normalised mutual information values (set to 1).
% 
% Alberto Liardi, 2024

function [U_X, U_Y, R, S, MI] = NMI_PID(Un_X, Un_Y, Red, Syn, MI)
    
    assert(length(MI)==length(Un_X) && length(MI)==length(Un_Y) ...
           && length(MI)==length(Red) && length(MI)==length(Syn))
    
    U_X = Un_X ./ MI;
    U_Y = Un_Y ./ MI;
    R = Red ./ MI;
    S = Syn ./ MI;
    MI(:) = 1;

end