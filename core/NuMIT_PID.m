%% NuMIT_PID
% 
% Returns the quantiles of a set of PID atoms with respect to a null model.
% 
% Description:
%   This function computes the quantile values for PID atoms with respect to a specified null model. 
%   The input PID atoms are compared against the correspondent null model and their quantile is computed.
% 
% Inputs:
%   UnX    - [L] array of Unique Information X atoms for which the quantiles are computed.
%   UnY    - [L] array of Unique Information Y atoms for which the quantiles are computed.
%   Red    - [L] array of Redundancy atoms for which the quantiles are computed.
%   Syn    - [L] array of Synergy atoms for which the quantiles are computed.
%   MI     - [L] array of total mutual information values.
%   model  - Struct with the following fields:
%            name    - "Gauss", "VAR", or "Discrete", name of the model for the null model computation.
%            S       - Number of sources.
%            T       - Number of targets.
%            n       - Number of iterations of the null model for each set of PID atoms.
%            p       - [L] array of model orders.
%            red_fun - Redundancy function to use: "MMI" (default), "DEP", or "CCS".
% 
% Outputs:
%   qUnX   - [L] array of quantile values for UnX atoms computed with the null model.
%   qUnY   - [L] array of quantile values for UnY atoms computed with the null model.
%   qRed   - [L] array of quantile values for Red atoms computed with the null model.
%   qSyn   - [L] array of quantile values for Syn atoms computed with the null model.
% 
% Alberto Liardi, 2024


function [qUnX, qUnY, qRed, qSyn] = NuMIT_PID(UnX, UnY, Red, Syn, MI, model)

    if isfield(model, 'p')
        % change length in size(,1)
        assert(length(MI)==length(model.p));
    else, model.p=ones(length(MI)); 
    end
    if ~isfield(model, 'n'), model.n=[]; end
    if ~isfield(model, 'S'), model.S=2; end
    if ~isfield(model, 'T'), model.T=1; end
    if ~isfield(model, 'red_fun'), model.red_fun="MMI"; end
    
    if model.name == "Discrete", r=1;
    else
        if model.red_fun=="MMI", r=1;
        elseif model.red_fun=="DEP", r=2;
        elseif model.red_fun=="CCS", r=3;
        end
    end
    
    nPIDs = cell(1,length(MI));
    if model.name == "Gauss"
        for l = 1:length(MI)
            null_pids = Null_model_Gauss(MI(l),model.S,model.T,model.n,model.red_fun); 
            if model.red_fun=="MMI"
                % NB: MMI assigns zero to either of the unique information, hence the null
                % distribution is the combination of the two (to avoid sequences of zeros)
                null_pids = [sum(null_pids(1:2,:,:),1); sum(null_pids(1:2,:,:),1); null_pids(3:4,:,:)];
            end
            nPIDs{l} = null_pids(:,:,r);
        end
    
    elseif model.name == "VAR"
        for l = 1:length(MI) 
            null_pids = Null_model_VAR(MI(l), model.S, model.p(l), model.n, ...
                                         model.red_fun, ["PID"]);
            if model.red_fun=="MMI"
                null_pids = [sum(null_pids(1:2,:,:),1); sum(null_pids(1:2,:,:),1); null_pids(3:4,:,:)];
            end
            nPIDs{l} = null_pids(:,:,r); 
        end

    elseif model.name == "Discrete"
        for l = 1:length(MI) 
            null_pids = Null_model_Discrete(MI(l),model.n);
            nPIDs{l} = [sum(null_pids(1:2,:),1); sum(null_pids(1:2,:),1); null_pids(3:4,:)]; 
        end

    else 
        error("Not a valid model type inserted, exiting."); 
    end
    
    qPIDs = zeros(4,length(MI));
    for m = 1:length(MI) 
        values = [UnX(m), UnY(m), Red(m), Syn(m)]';
        qPIDs(1:4,m) = CompQuantile(nPIDs{m},values);
    end
    
    qUnX = qPIDs(1,:);  qUnY = qPIDs(2,:);
    qRed = qPIDs(3,:);  qSyn = qPIDs(4,:);
       
end