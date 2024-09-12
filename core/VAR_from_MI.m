%% VAR_from_MI
% 
% Generates a random Vector Autoregressive (VAR) system with a specified total mutual information (MI).
% 
% Description:
%   This function generates a random VAR(p) model with a given number of sources and 
%   a specified total mutual information (MI) between the sources and targets. 
%   It optimises the spectral radius of the VAR coefficients to obtain the desired MI. 
% 
% Inputs:
%   MI      - Target total mutual information between sources and targets.
%   p       - Order of the VAR(p) model.
%   S       - Number of sources in the VAR system.
% 
% Outputs:
%   Coef    - Cell array of size [1, p] containing the coefficient matrices of the VAR system for each lag.
%   V       - [S, S] matrix representing the noise covariance (Wishart-distributed).
%   flag    - Scalar flag indicating the success (0) or failure (1) of the optimization process.
%   g       - Optimized spectral parameter for the system.
% 
% Notes:
%   - The function generates random initial VAR coefficients and a noise covariance matrix from a Wishart distribution.
%   - The optimization process adjusts the spectral radius of the VAR system to achieve the target MI using a Sigmoid transformation and the `fzero` function.
%   - If the optimization fails or returns NaN values, the function sets `flag` to 1 and returns. 
%   - The function includes a sub-function `fun` to compute the mutual information for a given spectral radius and another sub-function `Sigmoid` for the transformation.
%   - A warning message is printed if the optimization fails.
% 
% Alberto Liardi, 2024


function [Coef,V,flag,g] = VAR_from_MI(MI,p,S)
    
    opts = optimset('display', 'none');
    Coef = cell(1,p);
    % initial conditions for the optimisation
    x0 = [0 -1 2 -2 5 -5];
    flag=0;
    
    % random VAR matrix
    A = var_rand(S,p,rand(1));
    while any(abs(A(:))<1e-18)
        A = var_rand(S,p,rand(1));
    end
    
    % target conditional covariance
    % NB. for Wishart degrees of freedom must be > size - 1
    V = wishrnd(eye(S),S+1);
    
    % find g such that mutual information is MI_tot
    l = 1;
    %         t0 = tic();
    [g, ~, code, out] = fzero(@(x) fun(x,A,V,MI), x0(l), opts);
      
    while( code ~= 1 && l+1 <= length(x0) )
        l = l+1;
        [g, ~, code, out] = fzero(@(x) fun(x,A,V,MI), x0(l), opts);
    end
    if( code ~= 1 ), fprintf("optimisation failed "+sprintf('%01d',code)+"\n"); 
         flag=1; return; 
    end 
    
    % now do one more calculation
    gg = Sigmoid(g,1);
    [B,~] = specnorm(A,gg);
    for s = 1:p, Coef{s} = B(:,:,s); end

    if ~isfinite(fun(g,B,V,MI))
        disp("Optimisation function returned NaN value!");
        flag=1;
        return;
    end

end


function y = fun(x,A,V,I)

    x = Sigmoid(x,1);
    [BB,~] = specnorm(A,x);

    % sometimes for large p and specific (small) values of the spectral radius 
    % the Lyapunov equation does not have a unique solution
    try
        G = var_to_autocov(BB,V,0);
        MI_value = (0.5*log(det(G(:,:,1)))-0.5*log(det(V)))/log(2);
        if ~isreal(MI_value)
            % when the spectral radius is close to 1 the MI can become complex
            error("error in optimising the g!");
        end
        y = MI_value-I;

    catch 
        y=NaN;
        disp("NAN value encountered");

    end

end

function y = Sigmoid(x,M)
    
    y = M/(1+exp(-x));

end