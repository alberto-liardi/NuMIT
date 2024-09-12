%% PIDcheck
% 
% Validate the consistency of PID decomposition.
% 
% Description:
%   This function checks the consistency of PID values for each row of the input matrix. It verifies that
%   the both marginal and total mutual information are non-negative (threshold value is 1e-11).
% 
% Inputs:
%   PID    - [N,4] matrix where each row contains PID atoms in the order: 
%            [Unique X, Unique Y, Redundancy, Synergy].
% 
% Outputs:
%   None. The function throws an error if any of the checks fail.
% 
% Alberto Liardi, 2024

function [] = PIDcheck(PID)
    % values above -1e-11 are considered non-negative

    % in case of a matrix of PID atoms (each row consists of PID atoms for
    % a specific redundancy function)
    for s = 1:size(PID,1)

        Ux = PID(s,1);
        Uy = PID(s,2);
        R = PID(s,3);
        S = PID(s,4);
        
        MI = Ux + Uy + R + S;
        MIx = Ux + R;
        MIy = Uy + R;
        assert(MIx > -1e-11 && MIy > -1e-11, "Single mutual information is negative");
        assert(MI > max(MIx, MIy) - 5e-10, "MI is smaller than a singular MI");
%         assert(S > -1e-11, "Synergy is negative!");

    end
end