%% NuMIT_Phis
% 
% Returns the quantiles of Phi measures with respect to a null model.
% 
% Description:
%   This function computes the quantile values for Phi measures (Phi_WMS, Phi_R, Redundancy) 
%   with respect to a specified null model. It generates the null model based on a VAR model and computes the quantiles 
%   of the Phi measures by comparing them against those from the null model.
% 
% Inputs:
%   Phi_wms  - [L] array of Phi Whole-minus-sum (WMS) Phi measures for which the quantiles are computed.
%   Phi_R    - [L] array of Phi_r measures for which the quantiles are computed.
%   Phi_red  - [L] array of Redundancy measures for which the quantiles are computed.
%   MI       - [L] array of total mutual information values.
%   model    - Struct with the following fields:
%              name  - "VAR", name of the model (currently, only "VAR" is supported).
%              S     - Number of sources.
%              T     - Number of targets.
%              n     - Number of iterations of the null model for each set of Phi measures.
%              p     - [L] array of model orders (only one model order should be provided).
% 
% Outputs:
%   qPhi_wms - [L] array of quantile values for Phi_wms measures computed with the null model.
%   qPhi_R   - [L] array of quantile values for Phi_R measures computed with the null model.
%   qPhi_red - [L] array of quantile values for Phi_red measures computed with the null model.
% 
% Notes:
%   - VAR procedures only support model order p=1.
% 
% Alberto Liardi, 2024

function [qPhi_wms, qPhi_R, qPhi_red] = NuMIT_Phis(Phi_wms, Phi_R, Phi_red, MI, model)

if isfield(model, 'p')
    % change length in size(,1)
    assert(length(MI)==length(model.p));
else 
    model.p=ones(length(MI));
end
if ~isfield(model, 'n'), model.n=[]; end
if ~isfield(model, 'S'), model.S=2; end
if ~isfield(model, 'T'), model.T=1; end

nPhis = cell(1,length(MI));

if model.name == "VAR"
    for l = 1:length(MI)
        [~, nPhis{l}] = Null_model_VAR(MI(l),model.S,model.p(l),model.n,["Phis"]);
    end
else 
    error("Not a valid model type inserted, exiting."); 
end

qPhis = zeros(3,length(MI));
for m = 1:length(MI) 
    values = [Phi_wms(m), Phi_R(m), Phi_red(m)]';
    qPhis(1:3,m) = CompQuantile(nPhis{m},values);
end

qPhi_wms = qPhis(1,:);  qPhi_R = qPhis(2,:);  qPhi_red = qPhis(3,:);
