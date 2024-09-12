%% CompQuantile
% 
% Compute the quantiles of given values for input distributions.
% 
% Description:
%   This function computes the quantiles of the specified values for each
%   distribution provided in the input matrix. Each distribution is represented
%   by a row in the 'distr' matrix, and the quantile is computed for each value
%   in the 'value' matrix, which corresponds to the distributions in 'distr'.
% 
% Inputs:
%   distr     - [M,N] matrix where each row represents a distribution with N samples.
%   value     - [M,L] matrix where each row contains L values for which quantiles are computed.
% 
% Outputs:
%   quantile  - [M,L] matrix where each element represents the quantile of the corresponding
%               value in 'value' relative to the corresponding distribution in 'distr'.
% 
% Alberto Liardi, 2024


function quantile = CompQuantile(distr,value)

    assert(size(distr,1)==size(value,1), "Dimensions not consistent!");
    
    quantile = zeros(size(value));

    for l = 1:size(value,2)
        nless = sum(distr < value(:,l) - 1e-8,2);
        nequal = sum(abs(distr-value(:,l))<1e-8,2);
        quantile(:,l) = (nless + 0.5.*nequal) / length(rmmissing(distr(1,:)));
    end    
end
