%% EvolCoefFit
% 
% Estimate the coefficients of a Vector Autoregressive (VAR) model from time series data.
% 
% Description:
%   This function fits a VAR(p) model to the input time series data using a specified maximum model order.
%   It returns the estimated coefficients, model order, residual covariance, and residuals time series. 
% 
% Inputs:
%   X       - [channels, time points, trials] matrix of time series data.
%   momax   - Maximum model order for model order selection (optional, default is 15).
% 
% Outputs:
%   AA      - Cell array of size [1, morder] where each element is the [channels, channels] matrix of 
%             VAR coefficients for each lag.
%   morder  - Optimal VAR model p.
%   V       - [channels, channels] residual covariance matrix.
%   E       - [channels, time points, trials] matrix of residuals from the VAR model.
%
%
% The code was adapted from: https://github.com/SacklerCentre/MVGC2/blob/master/demo/mvgc_demo_var.m
% 
% Alberto Liardi, 2024

function [ AA,morder,V,E ] = EvolCoefFit(X, momax)

    % VAR model order estimation
    if ~exist('momax',     'var'), momax     = 15; end % maximum model order for model order selection
    if ~exist('seed',      'var'), seed      = 0;       end % random seed (0 for unseeded)  
    
    
    % Remove temporal mean and normalise by temporal variance.   
    X = demean(X,true);
    

    % Model order estimation:    
    % Calculate and plot VAR model order estimation criteria up to specified maximum model order.
%     ptic('\n*** tsdata_to_varmo... ');
    [moaic,mobic,mohqc,molrt] = tsdata_to_varmo(X,momax,'LWR',[],[],[]);
%     ptoc;
    
    % Select and report VAR model order.
    morder = moselect(sprintf('VAR model order selection (max = %d)',momax), 'LRT','AIC',moaic,'BIC',mobic,'HQC',mohqc,'LRT',molrt);
    assert(morder > 0,'selected zero model order! GCs will all be zero!');
%     if morder >= momax, fprintf(2,'*** WARNING: selected maximum model order (may have been set too low)\n'); end
    

    % VAR model estimation:
    % Estimate VAR model of selected order from data.
%     ptic('\n*** tsdata_to_var... ');
    [A,V,E] = tsdata_to_var(X,morder,'LWR');
%     ptoc;
    
    % Check for failed regression
    assert(~isbad(A),'VAR estimation failed - bailing out');
    info = var_info(A,V,0);
    assert(~info.error,'VAR error(s) found - bailing out');
    
    % saving coefficient matrices in cells for ease of use
    AA = cell(1,morder);
    for l = 1:morder
        AA{l} = A(:,:,l);
    end

end