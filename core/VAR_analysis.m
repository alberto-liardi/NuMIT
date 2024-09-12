%% VAR_analysis
% 
% Fits time series data to a VAR(p) model and computes PID, PhiID, and Phi measures.
% 
% Description:
%   This function fits time series data to a Vector Autoregressive (VAR) model of order p 
%   and computes information measures such as Partial Information Decomposition (PID), 
%   Integrated Information Decomposition (PhiID), and Integrated Information (Phi) metrics. 
%   The analysis is run over multiple iterations with random selections of channels and epochs. 
% 
% Inputs:
%   data           - {E} cell array of [T,C] matrices. Each cell represents an epoch containing 
%                    time series data with C channels and T time points.
%   subject_info   - Struct with metadata about the subject, including:
%                    * fs        - Sampling frequency of the data (required).
%                    * condition - Subject condition (e.g., placebo, awake).
%   par            - Struct with the analysis parameters, including:
%                    * channels  - Number of random channels for the VAR model (default: 2).
%                    * epochs    - Number of random epochs for the VAR model (default: 30).
%                    * runs      - Number of iterations of the procedure (default: 100).
%                    * mmorder   - Maximum model order for the VAR model (default: 30).
%                                  !! For PhiID and Phis, mmorder is overwritten and set to 1.
%                                  !! If p > 1, only "PID" can be used as a method.
%                    * red_fun   - Redundancy functions to use: "MMI" (default), "DEP", "CCS", or "all" (for all functions).
%                                  !! Only useful for PID; for PhiID and Phis, only MMI is used.
%   methods        - Array containing the information metrics to compute. Available options: 
%                    "PID", "PhiID", "Phis" (default: ["PID", "PhiID", "Phis"]).
%                    !! Choose ["PID"] if p > 1.
% 
% Outputs:
%   Information    - Struct containing the computed information measures:
%                    * par          - The analysis parameters (see Inputs).
%                    * subject_info - Metadata of the subject (see Inputs).
%                    * PIDs         - (If "PID" is selected) Cell array of length 'runs' containing [3,4] matrices 
%                                     for each run. Rows correspond to MMI-, DEP-, and CCS-PID atoms: 
%                                     [UniqueX, UniqueY, Redundancy, Synergy].
%                    * PhiIDs       - (If "PhiID" is selected) Cell array of length 'runs' containing structs 
%                                     with PhiID atoms for each run.
%                    * Phi_WMS      - (If "Phis" is selected) Array of length 'runs' containing 
%                                     Phi whole-minus-sum values.
%                    * Phi_R        - (If "Phis" is selected) Array of length 'runs' containing Phi_R values.
% 
% Notes:
%   - The time series data is filtered and downsampled to 100 Hz before fitting the VAR model.
%   - If PhiID or Phis is selected as a method, the model order is forced to p = 1.
%   - The function fits a VAR model and computes the information measures for multiple randomly selected channels and epochs.
% 
% Alberto Liardi, 2024
                

function Information = VAR_analysis(data, subject_info, par, methods)

    % setting parameters of the data
    tot_epochs = length(data);
    tot_channels = size(data{1},2);
    tot_period = size(data{1},1);
    
    % setting parameter for the analyses
    if ~isfield(par, 'channels'), par.channels = 2; end
    if ~isfield(par, 'epochs'), par.epochs = 30; end, if tot_epochs<par.epochs, par.epochs=tot_epochs; end
    if ~isfield(par, 'runs'), par.runs = 100; end
    if ~isfield(par, 'mmorder'), par.mmorder = 30; end
    if ~exist('methods','var'), methods = ["PID"]; end

    if ismember("PID", methods)
        % choosing redundancy function (MMI/DEP/CCS)
        if ~isfield(par, 'red_fun'), par.red_fun = ["MMI"]; end
        PIDs = cell(1, par.runs);
    end

    if ismember("PhiID", methods)
        par.mmorder = 1;
        PhiIDs = cell(1, par.runs);
    end
    
    if ismember("Phis", methods)
        par.mmorder = 1;
        Phi_WMS = zeros(1,par.runs);
        Phi_R = zeros(1,par.runs);
    end
    
    % repeating the procedure par.runs times
    for N = 1:1:par.runs
    
        % choosing 'par.channels' random channels and 'par.epochs' random epochs
        c = randperm(tot_channels, par.channels);
        e = randperm(tot_epochs, par.epochs);   
    
        % assembling the time series
        Z = zeros(par.channels, tot_period, par.epochs);
        for j = 1:par.channels
            for n = 1:par.epochs
                Z(j,:,n) = data{e(n)}(:,c(j))';
            end
        end
    
        % filtering and downsampling to 100Hz
        down_s = floor(subject_info.fs/200);
        Y = zeros(par.channels, tot_period/down_s, par.epochs);
        for j = 1:par.channels
            for n = 1:par.epochs
                Y(j,:,n) = LowPassFilter(Z(j,:,n), down_s, subject_info.fs);
            end
        end
            

        %%%%%%% VAR fitness %%%%%%%
        % fitting the evolution coefficient and (optional) residuals
        % NB: the model order is fixed at p=1
        [A,p,V,~] = EvolCoefFit(Y, par.mmorder);
                

        %%%%%%% computing information measures %%%%%%%     
        
        % computing PID atoms   
        if ismember("PID", methods)
            [PIDs{N}, Gammas] = PID_VAR_calculator(p,A,V,par.channels/2, ...
                                    par.channels/2, par.red_fun);
        end

        % computing PhiID atoms
        if ismember("PhiID", methods)
            [PhiIDs{N}, Gammas] = PhiID_VAR_calculator(p,A,V,par.channels/2,par.channels/2);
        end
        
        % computing Phi_WMS and Phi_R
        if ismember("Phis", methods)
            % if gamma haven't been computed already, do it now
            if ~exist('Gammas','var'), Gammas = get_Gamma(A,V,p); end
            MI = (0.5*log(det(Gammas{1}))-0.5*log(det(V)))/log(2);
            [Phi_WMS(N), Phi_R(N)] = Phi_VAR_calculator(MI,Gammas);
        end
    end
    
    % save everything in the Information struct

    Information.subject_info = subject_info;
    Information.parameters = par;
    
    if ismember("PID", methods)
        Information.PIDs = PIDs;
    end

    if ismember("PhiID", methods)
        Information.PhiIDs = PhiIDs;
    end

    if ismember("Phis", methods)
        Information.phi_WMS = Phi_WMS;
        Information.phi_R = Phi_R;
    end
end


% return gamma autocovariances from VAR parameters in cell array 
function Gamma = get_Gamma(A,V,p)
    G = var_to_autocov(cat(3, A{:}),V,p);
    Gamma = cell(1,p+1);
    for s = 1:p+1
        Gamma{s} = G(:,:,s);
    end
end