%% DiscreteMarginalMI
% 
% Compute the marginal mutual information for a discrete system with a specified logic gate.
% 
% Description:
%   This function calculates the marginal mutual information (MI) of a discrete system where the target
%   is given by T = f(X,Y) + pe. The function computes the mutual information based on the type of
%   2-bit logic gate specified. The gates supported are XOR, XOR2, XOR3, OR, OR2, OR3, and OR4. The probabilities
%   provided describe the joint distribution of X and Y, and `pe` is the probability of flipping the outcome of
%   f(X,Y). `HT` is the entropy of the target distribution.
% 
% Inputs:
%   ps   - 1x4 array of probabilities: [P(x=0,y=0), P(x=0,y=1), P(x=1,y=0), P(x=1,y=1)].
%   pe   - Probability of flipping the outcome of f(X,Y) (from 0 to 1 and vice versa).
%   HT   - Entropy of the target.
%   gate - Type of 2-bit logic gate to use ('XOR', 'XOR2', 'XOR3', 'OR', 'OR2', 'OR3', 'OR4').
% 
% Outputs:
%   MIx  - Marginal mutual information with respect to variable X.
%   MIy  - Marginal mutual information with respect to variable Y.
% 
% Alberto Liardi, 2024

function [MIx, MIy] = DiscreteMarginalMI(ps,pe,HT,gate)

    if gate=="XOR"
        [MIx, MIy] = XORgate(ps,pe,HT);
    elseif gate=="XOR2"
        [MIx, MIy] = XOR2gate(ps,pe,HT);
    elseif gate=="XOR3"
        [MIx, MIy] = XOR3gate(ps,pe,HT);
    elseif gate=="OR"
        [MIx, MIy] = ORgate(ps,pe,HT);
    elseif gate=="OR2"
        [MIx, MIy] = OR2gate(ps,pe,HT);
    elseif gate=="OR3"
        [MIx, MIy] = OR3gate(ps,pe,HT);
    elseif gate=="OR4"
        [MIx, MIy] = OR4gate(ps,pe,HT);
    end
end
    

%%%%% XOR gates
    
function [MIx, MIy] = XORgate(ps,pe,HT)

    % unpack probabilities
    cps = num2cell(ps);
    [p00,p01,p10,p11] = cps{:};

    % single probabilities
    px0 = p00 + p01;
    px1 = 1-px0;
    
    py0 = p00 + p10;
    py1 = 1-py0;

    % conditional probabilities on Z for X
    pZ0_x0 = p00/px0;
    pZ1_x0 = p01/px0;
    pZ0_x1 = p11/px1;
    pZ1_x1 = p10/px1;
    
    % conditional probabilities on Z for Y
    pZ0_y0 = p00/py0;
    pZ1_y0 = p10/py0;
    pZ0_y1 = p11/py1;
    pZ1_y1 = p01/py1;

    [MIx, MIy] = computeMarginals(pZ0_x0, pZ1_x0, pZ0_x1, pZ1_x1, ...
                                  pZ0_y0, pZ1_y0, pZ0_y1, pZ1_y1, ...
                                  px0, px1, py0, py1, pe, HT);
    
end

function [MIx, MIy] = XOR2gate(ps,pe,HT)

    % unpack probabilities
    cps = num2cell(ps);
    [p00,p01,p10,p11] = cps{:};

    % single probabilities
    px0 = p00 + p01;
    px1 = 1-px0;
    
    py0 = p00 + p10;
    py1 = 1-py0;

    % conditional probabilities on Z for X
    pZ0_x0 = 1;
    pZ1_x0 = 0;
    pZ0_x1 = 0;
    pZ1_x1 = 1;
    
    % conditional probabilities on Z for Y
    pZ0_y0 = p00/py0;
    pZ1_y0 = p10/py0;
    pZ0_y1 = p01/py1;
    pZ1_y1 = p11/py1;

    [MIx, MIy] = computeMarginals(pZ0_x0, pZ1_x0, pZ0_x1, pZ1_x1, ...
                                  pZ0_y0, pZ1_y0, pZ0_y1, pZ1_y1, ...
                                  px0, px1, py0, py1, pe, HT);
    
end

function [MIx, MIy] = XOR3gate(ps,pe,HT)

    % unpack probabilities
    cps = num2cell(ps);
    [p00,p01,p10,p11] = cps{:};

    % single probabilities
    px0 = p00 + p01;
    px1 = 1-px0;
    
    py0 = p00 + p10;
    py1 = 1-py0;

    % conditional probabilities on Z for X
    pZ0_x0 = p00/px0;
    pZ1_x0 = p01/px0;
    pZ0_x1 = p10/px1;
    pZ1_x1 = p11/px1;
    
    % conditional probabilities on Z for Y
    pZ0_y0 = 1;
    pZ1_y0 = 0;
    pZ0_y1 = 0;
    pZ1_y1 = 1;

    [MIx, MIy] = computeMarginals(pZ0_x0, pZ1_x0, pZ0_x1, pZ1_x1, ...
                                  pZ0_y0, pZ1_y0, pZ0_y1, pZ1_y1, ...
                                  px0, px1, py0, py1, pe, HT);
    
end

%%%%% OR gates now

function [MIx, MIy] = ORgate(ps,pe,HT)

    % unpack probabilities
    cps = num2cell(ps);
    [p00,p01,p10,p11] = cps{:};

    % single probabilities
    px0 = p00 + p01;
    px1 = 1-px0;
    
    py0 = p00 + p10;
    py1 = 1-py0;

    % X
    % conditional probabilities on Z
    pZ0_x0 = p00/px0;
    pZ1_x0 = p01/px0;
    pZ0_x1 = 0;
    pZ1_x1 = 1;
    
    % Y
    % conditional probabilities on Z
    pZ0_y0 = p00/py0;
    pZ1_y0 = p10/py0;
    pZ0_y1 = 0;
    pZ1_y1 = 1;

    [MIx, MIy] = computeMarginals(pZ0_x0, pZ1_x0, pZ0_x1, pZ1_x1, ...
                                  pZ0_y0, pZ1_y0, pZ0_y1, pZ1_y1, ...
                                  px0, px1, py0, py1, pe, HT);

end

function [MIx, MIy] = OR2gate(ps,pe,HT)

    % unpack probabilities
    cps = num2cell(ps);
    [p00,p01,p10,p11] = cps{:};

    % single probabilities
    px0 = p00 + p01;
    px1 = 1-px0;
    
    py0 = p00 + p10;
    py1 = 1-py0;

    % X
    % conditional probabilities on Z
    pZ0_x0 = p01/px0;
    pZ1_x0 = p00/px0;
    pZ0_x1 = 0;
    pZ1_x1 = 1;
    
    % Y
    % conditional probabilities on Z
    pZ0_y0 = 0;
    pZ1_y0 = 1;
    pZ0_y1 = p01/py1;
    pZ1_y1 = p11/py1;

    [MIx, MIy] = computeMarginals(pZ0_x0, pZ1_x0, pZ0_x1, pZ1_x1, ...
                                  pZ0_y0, pZ1_y0, pZ0_y1, pZ1_y1, ...
                                  px0, px1, py0, py1, pe, HT);

end

function [MIx, MIy] = OR3gate(ps,pe,HT)

    % unpack probabilities
    cps = num2cell(ps);
    [p00,p01,p10,p11] = cps{:};

    % single probabilities
    px0 = p00 + p01;
    px1 = 1-px0;
    
    py0 = p00 + p10;
    py1 = 1-py0;

    % X
    % conditional probabilities on Z
    pZ0_x0 = 0;
    pZ1_x0 = 1;
    pZ0_x1 = p10/px1;
    pZ1_x1 = p11/px1;
    
    % Y
    % conditional probabilities on Z
    pZ0_y0 = p10/py0;
    pZ1_y0 = p00/py0;
    pZ0_y1 = 0;
    pZ1_y1 = 1;

    [MIx, MIy] = computeMarginals(pZ0_x0, pZ1_x0, pZ0_x1, pZ1_x1, ...
                                  pZ0_y0, pZ1_y0, pZ0_y1, pZ1_y1, ...
                                  px0, px1, py0, py1, pe, HT);

end

function [MIx, MIy] = OR4gate(ps,pe,HT)

    % unpack probabilities
    cps = num2cell(ps);
    [p00,p01,p10,p11] = cps{:};

    % single probabilities
    px0 = p00 + p01;
    px1 = 1-px0;
    
    py0 = p00 + p10;
    py1 = 1-py0;

    % X
    % conditional probabilities on Z
    pZ0_x0 = 0;
    pZ1_x0 = 1;
    pZ0_x1 = p11/px1;
    pZ1_x1 = p10/px1;
    
    % Y
    % conditional probabilities on Z
    pZ0_y0 = 0;
    pZ1_y0 = 1;
    pZ0_y1 = p11/py1;
    pZ1_y1 = p01/py1;

    [MIx, MIy] = computeMarginals(pZ0_x0, pZ1_x0, pZ0_x1, pZ1_x1, ...
                                  pZ0_y0, pZ1_y0, pZ0_y1, pZ1_y1, ...
                                  px0, px1, py0, py1, pe, HT);

end


%%%%% Marginal MI calcuator

function [MIx, MIy] = computeMarginals(pZ0_x0, pZ1_x0, pZ0_x1, pZ1_x1, ...
                                       pZ0_y0, pZ1_y0, pZ0_y1, pZ1_y1, ...
                                       px0, px1, py0, py1, pe, HT)

    % X
    % conditional probabilities on T
    pT0_x0 = pZ0_x0*(1-pe) + pZ1_x0*pe;
    pT1_x0 = pZ1_x0*(1-pe) + pZ0_x0*pe; % 1-pT0_x0;
    pT0_x1 = pZ0_x1*(1-pe) + pZ1_x1*pe;
    pT1_x1 = pZ1_x1*(1-pe) + pZ0_x1*pe; % 1-pT0_x1;
    
    % conditional entropies
    HT_x0 = -px0*(pT0_x0*log(pT0_x0) + pT1_x0*log(pT1_x0));
    HT_x1 = -px1*(pT0_x1*log(pT0_x1) + pT1_x1*log(pT1_x1));
    HTcondX = HT_x0 + HT_x1;

    % Y
    % conditional probabilities on T
    pT0_y0 = pZ0_y0*(1-pe) + pZ1_y0*pe;
    pT1_y0 = pZ1_y0*(1-pe) + pZ0_y0*pe; % 1-pT0_y0;
    pT0_y1 = pZ0_y1*(1-pe) + pZ1_y1*pe;
    pT1_y1 = pZ1_y1*(1-pe) + pZ0_y1*pe; % 1-pT0_y1;
    
    % conditional entropies
    HT_y0 = -py0*(pT0_y0*log(pT0_y0) + pT1_y0*log(pT1_y0));
    HT_y1 = -py1*(pT0_y1*log(pT0_y1) + pT1_y1*log(pT1_y1));
    HTcondY = HT_y0 + HT_y1;
    
    % mutual information
    MIx = HT - HTcondX;
    MIy = HT - HTcondY;

    if MIx<0 || MIy<0
        error("Marginal MI are negative!");
    end

end