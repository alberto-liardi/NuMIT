%% Discrete_from_MI
% 
% Generates a discrete system of the form T = Z + epsilon for a specified mutual information (MI) value.
% 
% Description:
%   This function computes the value of the discrete noise parameter epsilon such that the mutual 
%   information between the system components matches the target MI value. 
%   The system is defined by T = Z + epsilon, where Z = gate(X,Y) is a Boolean function of X and Y. 
%   The function returns the noise parameter epsilon and a flag indicating whether the optimization succeeded.
% 
% Inputs:
%   MI      - Target mutual information value for the discrete system.
%   ps      - [1x4] vector of joint probabilities of X and Y.
%   boolfun - Function handle for the Boolean gate, which computes the probability of Z given ps.
% 
% Outputs:
%   pe      - Discrete noise parameter epsilon that matches the specified MI.
%   flag    - Flag indicating whether the optimization succeeded (0) or failed (1).
% 
% 
% Alberto Liardi, 2024


function [pe, flag] = Discrete_from_MI(MI,ps,boolfun)

    assert(abs(sum(ps)-1)<1e-6, "Probabilities do not sum up to 1!");
    flag=0;

    % set parameters for the optimisation
    x0 = [1 -1 2 -2 5 -5];
    opts = optimset('display', 'none');

    % get the probability of Z=gate(X,Y)
    [pz0,pz1] = boolfun(ps);

    % for some parameters of the Dirichlet the MI value may not be reached,
    % check with the maximum of MI (at x=1 or x=0) and skip the step if needed
    if -pz1*log(pz1)-pz0*log(pz0)<MI, disp("Sampled parameters not valid, skipping step!"); 
        flag=1; 
        pe=NaN; 
        return; 
    end
    
    % find the discrete noise epsilon to get the value of MI
    l = 1;
    [f, ~, code, ~] = fzero(@(x) discrfun(x, pz0, pz1, MI), x0(l), opts);
    
    while( code ~= 1 && l+1 <= length(x0) )
        l = l+1;
        [f, ~, code, ~] = fzero(@(x) discrfun(x, pz0, pz1, MI), x0(l), opts);
    end
    if( code ~= 1 )
        flag = 1; 
        fprintf("optimisation failed!\n");
        pe=NaN;
        return; 
    end 
    
    % now that the optimisation is finished, do one more calculation
    pe = Sigmoid(f,0.5);

    if ~isfinite(discrfun(f, pz0, pz1, MI))
        flag = 1;
        disp("Optimisation function returned NaN value!");
    end

end


function y = discrfun(x, pz0, pz1, MI_value)
    
    x = Sigmoid(x,0.5);

    % can calculate the entropy of the target T=Z+epsilon manually:
    pt0 = pz0*(1-x) + pz1*x;
    pt1 = pz1*(1-x) + pz0*x; % 1-pt0;
    
    HT = -pt0*log(pt0)-pt1*log(pt1);
    HTcondZ = -x*log(x)-(1-x)*log(1-x);
    I = HT - HTcondZ;
    
    y = MI_value - I;

end

function y = Sigmoid(x,M)
    
    y = M/(1+exp(-x));

end