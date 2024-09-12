%% PID_MMI_Discrete
%
% Computes the MMI-based Partial Information Decomposition on a discrete system.
%
% Description:
%   This function calculates the PID atoms (Unique X, Unique Y, Redundancy, Synergy) for a discrete system
%   of the form T = Z + epsilon, where Z is the output of a specified 2-input logic gate applied to X and Y,
%   and epsilon is a binary noise that flips Z with probability pe. The function returns the PID atoms and the 
%   total mutual information of the system.
%
% Inputs:
%   ps    - [4x1] array of representing the probabilities [p(X=0, Y=0); p(X=0, Y=1); p(X=1, Y=0); p(X=1, Y=1)].
%   pe    - Probability of flipping the output Z.
%   gate  - String specifying the 2-input logic gate used in the system. 
%           Options include "XOR", "XOR2", "XOR3", "OR", "OR2", "OR3", "OR4".
%
% Outputs:
%   UnX   - Unique information from X.
%   UnY   - Unique information from Y.
%   Red   - Redundant information shared between X and Y.
%   Syn   - Synergistic information between X and Y.
%   MI    - Total mutual information of the system.
%
% Alberto Liardi, 2024


function [UnX, UnY, Red, Syn, MI] = PID_MMI_Discrete(ps, pe, gate)

    if gate == "XOR"
        boolfun = @XORgate;
    elseif gate == "XOR2"
        boolfun = @XOR2gate;
    elseif gate == "XOR3"
        boolfun = @XOR3gate;
    elseif gate == "OR"
        boolfun = @ORgate;
    elseif gate == "OR2"
        boolfun = @OR2gate;
    elseif gate == "OR3"
        boolfun = @OR3gate;
    elseif gate == "OR4"
        boolfun = @OR4gate;
    end 
    
    [pz0,pz1] = boolfun(ps);
    
    pt0 = pz0*(1-pe) + pz1*pe;
    pt1 = pz1*(1-pe) + pz0*pe; % 1-pt0;
    
    HT = -pt0*log(pt0)-pt1*log(pt1);
    HTcondZ = -pe*log(pe)-(1-pe)*log(1-pe);
    MI = HT - HTcondZ;
    
    [MIx, MIy] = DiscreteMarginalMI(ps,pe,HT,gate);
    
    % now calculate PID atoms (MMI definition)
    Red = min(MIx, MIy);
    UnX = MIx - Red;
    UnY = MIy - Red;
    Syn = MI - UnX - UnY - Red;

end


function [pz0,pz1] = XORgate(ps)
    pz0 = ps(1)+ps(4);
    pz1 = 1-pz0;
end

function [pz0,pz1] = XOR2gate(ps)
    pz0 = ps(1)+ps(2);
    pz1 = 1-pz0;
end

function [pz0,pz1] = XOR3gate(ps)
    pz0 = ps(1)+ps(3);
    pz1 = 1-pz0;
end

function [pz0,pz1] = ORgate(ps)
    pz0 = ps(1);
    pz1 = 1-ps(1);
end

function [pz0,pz1] = OR2gate(ps)
    pz0 = ps(2);
    pz1 = 1-pz0;
end

function [pz0,pz1] = OR3gate(ps)
    pz0 = ps(3);
    pz1 = 1-pz0;
end

function [pz0,pz1] = OR4gate(ps)
    pz0 = ps(4);
    pz1 = 1-pz0;
end