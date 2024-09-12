%% PID_MMI_Gaussian
% 
% Compute Partial Information Decomposition (PID) components for a Gaussian system using
% the Minimal Mutual Information (MMI) redundancy definition.
% 
% Description:
%   This function calculates the MMI-PID atoms between two sources and a target 
%   in a Gaussian system. Sources and target variables can be multivariate.
% 
% Inputs:
%   Sigma  - [S+T, S+T] covariance matrix of the system.
%   S      - Number of sources in the system.
%   T      - Number of targets in the system.
% 
% Outputs:
%   UnX    - Unique information from the first source.
%   UnY    - Unique information from the second source.
%   Red    - Redundancy.
%   Syn    - Synergy.
%   MI_S1  - Marginal mutual information between the first source and the target.
%   MI_S2  - Marginal mutual information between the second source and the target.
%   MI     - Total mutual information between the sources and the target.
% 
% Alberto Liardi, 2024


function [UnX, UnY, Red, Syn, MI_S1, MI_S2, MI] = PID_MMI_Gaussian(Sigma,S,T)

    [~,pd] = chol(Sigma);
    assert(pd==0, "The covariance matrix provided is not positive definite!");
    assert(size(Sigma,1)==S+T, "Mismatch between number of sources/targets and covariance dimensions.");

    Sigma_t = Sigma(end-T+1:end,end-T+1:end);
    Sigma_s = Sigma(1:S,1:S);
    Sigma_cross = Sigma(1:S,S+1:end)';

    MI = entropy(Sigma_t)+entropy(Sigma_s)-entropy(Sigma);
            
    % S1
    Sigma_TS1 = [Sigma_t, Sigma_cross(:,1:S/2);
                 Sigma_cross(:,1:S/2)', Sigma_s(1:S/2,1:S/2)];
    MI_S1 = entropy(Sigma_t)+entropy(Sigma_s(1:S/2,1:S/2)) ...
            - entropy(Sigma_TS1);
    
    % S2
    Sigma_TS2 = [Sigma_t, Sigma_cross(:,S/2+1:end);
                 Sigma_cross(:,S/2+1:end)', Sigma_s(S/2+1:end,S/2+1:end)];
    MI_S2 = entropy(Sigma_t)+entropy(Sigma_s(S/2+1:end,S/2+1:end)) ...
            - entropy(Sigma_TS2);
    
    
    % now calculate PID atoms with MMI definition
    Red = min(MI_S1, MI_S2);
    UnX = MI_S1 - Red;
    UnY = MI_S2 - Red;
    Syn = MI - UnY - UnX - Red;

end


function H = entropy(Sig)
    
    arg = det(2*pi*exp(1)*Sig);    
    H = 0.5*log(arg);

end