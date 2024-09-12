%% Gaussian_from_MI
% 
% Generate a random Gaussian system with a specified total Mutual Information (MI) between sources and targets.
% 
% Description:
%   This function generates a covariance matrix for a Gaussian system with a specified mutual information
%   between the sources and targets. The function accepts source and target covariance matrices or the linear coefficients. 
%   If these are not provided, they are randomly sampled from a Wishart distribution. 
%   The function returns the covariance matrix of the system, the optimisation status, 
%   the linear coefficients, and a scaling factor that controls the mutual information.
% 
% Inputs:
%   MI        - Scalar value specifying the desired total mutual information between the sources and targets.
%   S         - Number of sources in the system (default is 2).
%   T         - Number of targets in the system (default is 1).
%   A         - [T,S] matrix of linear coefficients relating sources to targets (optional).
%   Sigma_s   - [S,S] covariance matrix of the sources (optional).
%   Sigma_u   - [T,T] conditional covariance matrix of the targets given the sources (optional).
% 
% Outputs:
%   Sigma     - [S+T, S+T] covariance matrix of the full system.
%   flag      - Scalar flag indicating success (0) or failure (1) of the mutual information optimization.
%   A         - [T,S] matrix of linear coefficients used in the system.
%   g         - Scalar scaling factor applied to the conditional covariance matrix.
% 
% Alberto Liardi, 2024


function [Sigma,flag,A,g] = Gaussian_from_MI(MI,S,T,A,Sigma_s,Sigma_u)

    if (~exist('S','var') || isempty(S)), S = 2; end
    if ~exist('T','var') || isempty(T), T = 1; end

    x0 = [1 2 0 5 10 20 50];
    opts = optimset('display', 'none');        
    flag=0;
    
    % sample the matrices
    % NB: Wishart degrees of freedom must be > size - 1
    
    if ~exist('Sigma_s','var') || isempty(Sigma_s)
        % sample source covariance
        Sigma_s = wishrnd(eye(S),S+1);
    end
    if ~exist('A','var') || isempty(A)
        % sample linear coefficients
        A = normrnd(0,1,T,S);
    end
    if ~exist('Sigma_u','var') || isempty(Sigma_u)
        % sample target conditional covariance
        Sigma_u = wishrnd(eye(T),T+1);
    end
    
    % find g such that mutual information is MI_tot
    l = 1;
    [alpha, ~, code] = fzero(@(x) fun(x,Sigma_u,A,Sigma_s,MI), x0(l), opts);
    
    while( (code == -3 || alpha<0) && l+1 <= length(x0) )
        l = l+1;
        [alpha, ~, code] = fzero(@(x) fun(x,Sigma_u,A,Sigma_s,MI), x0(l), opts);
    end
    if( code == -3 || alpha<0 ), disp("problem with MI"); flag=1; return; end 
    
    g = alpha^(-1);
    
    % check
    % 0.5*log(det(2*pi*exp(1)*(A*Sigma_s*A'+g*Sigma_u))) ...
    % - 0.5*log(det(2*pi*exp(1)*(g*Sigma_u)))
    
    % get the covariance matrix
    Sigma_t = A * Sigma_s * A.' + g*Sigma_u;
    
    % calculate cross covariance terms
    Sigma_cross = A*Sigma_s;

    Sigma = [Sigma_s, Sigma_cross';
                  Sigma_cross, Sigma_t];

    % check the result of the optimisation
    I_check = entropy(Sigma_t)+entropy(Sigma_s)-entropy(Sigma);

    if(abs(MI-I_check)>1e-6)
        flag=1; 
        disp("Error in the optimisation of the MI.");
%         disp(abs(MI-I_check)); error("problem with MI");
    end

end


function H = entropy(Sigma)
    
    arg = det(2*pi*exp(1)*Sigma);    
    H = 0.5*log(arg);

end

function y = fun(x,Sigma_u,A,Sigma_s,I)

    T = length(Sigma_u);
    y = det(eye(T)+ x * (Sigma_u \ A * Sigma_s' * A.'))-exp(2*I);

end