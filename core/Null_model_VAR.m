%% Null_model_VAR
% 
% Generates the null distribution of PID, PhiID, and Phis measures for VAR(p) systems with fixed total mutual information (MI).
%
% Description:
%   This function generates the null distribution of PID atoms, PhiID atoms, and integrated information measures
%   (Phi_WMS, Phi_R, and red_phi) for a VAR(p) model with a specified total mutual information (MI). The function 
%   computes the desired metrics based on the provided 'methods' flags.
%
% Inputs:
%   MI         - Target mutual information value for the VAR system.
%   S          - Number of sources in the system (must be even).
%   p          - Order of the VAR model (default is 1). PhiID and Phis only support p=1.
%   N_runs     - Number of runs for the simulation (default is 1,000).
%   red_fun    - Redundancy functions to use: "MMI" (default), "DEP", "CCS", or "all" (for all functions).
%                Only useful for PID; PhiID and Phis only support MMI.
%   methods    - Array containing the information metrics to compute: "PID", "PhiID", "Phis" (default: ["PID"]).
%                Choose ["PID"] for a VAR(p) model with p>1.
%
% Outputs:
%   PIDs       - [4, N_runs, 3] matrix where each slice represents PID atoms computed using different
%                redundancy functions:
%                [:,:,1] MMI definition
%                [:,:,2] DEP definition
%                [:,:,3] CCS definition
%                Each of the N_runs columns contains PID atoms in the order:
%                [Unique X; Unique Y; Redundancy; Synergy]
%                If redundancies are not selected, the matrix will contain zeros.
%   PhiIDs     - Cell array of size [1, N_runs], where each cell contains the PhiID atoms for the corresponding run.
%   Phis       - [3, N_runs] matrix where each row contains integrated information measures:
%                [Phi_WMS; Phi_R; red_phi].
%                phi_WMS - phi Whole Minus-Sum
%                phi_R   - phi_R (phi_WMS+red_phi)
%                red_phi  - Redundancy (with MMI definition)
%
% Notes:
%   - The function verifies the results and skips runs if the MI value is not reached or if errors occur during calculations.
%   - The PhiID and Phis measures only support VAR models of order p=1.
%   - The output matrix 'PIDs' has three slices for PID atoms, even if only one redundancy function is selected. 
%     In such case all the others will be zero.
% 
% Alberto Liardi, 2024

function [PIDs, PhiIDs, Phis] = Null_model_VAR(MI,S,p,N_runs,red_fun,methods)
    
    if ~exist('S','var') || isempty(S), S = 2; end
    if ~exist('p','var') || isempty(p), p = 1; end
    if ~exist('N_runs','var') || isempty(N_runs), N_runs = 1e4; end
    if ~exist('red_fun','var') || isempty(red_fun), red_fun = ["MMI"]; end
    if ~exist('methods','var'), methods = ["PID"]; end

    if mod(S,2)~=0, error("Provide even number of sources!"); end

    if ismember("PID", methods)
        rinds=[];
        if ismember("MMI",red_fun), rinds(end+1)=1; end
        if ismember("DEP",red_fun), rinds(end+1)=2; end
        if ismember("CCS",red_fun), rinds(end+1)=3; end
        if ismember("all",red_fun), rinds=1:3; end
    end

    if ismember("PhiID", methods) || ismember("Phis", methods)
        assert(p==1, "Invalid model order for PhiID or Phis, exiting.");
    end

    PIDs = NaN(4,N_runs,3);
    PhiIDs = cell(1, N_runs);
    phi_WMS = NaN(1,N_runs); phi_R = NaN(1,N_runs); red_phi = NaN(1,N_runs);

    nerr=0;
    
    for j = 1:N_runs

%         if mod(j,N_runs/10)==0, disp(j); end

        % generate random VAR model with specified MI
        [A,V,check] = VAR_from_MI(MI,p,S);
        
        if(check==1)
            disp("problem with MI optimisation"); 
            nerr = nerr+1;
            continue
        end

        if ismember("PID", methods)
            PID = PID_VAR_calculator(p,A,V,S/2,S/2,red_fun);
            
            % check MI optimisation
            for r = rinds
                if check_MI_PID(PID(r,:), MI)==0, nerr=nerr+1; continue; end
            end
                    
            PIDs(:,j,:) = PID';

        end

        if ismember("PhiID", methods)
            PhiID = PhiID_VAR_calculator(p,A,V,S/2,S/2);
            if check_MI_PhiID(PhiID, MI)==0, nerr=nerr+1; continue; end
            PhiIDs{j} = PhiID;
        end
        
        if ismember("Phis", methods)
            Gs = get_Gamma(A,V,p);
            MI_calc = (0.5*log(det(Gs{1}))-0.5*log(det(V)))/log(2);
            [phi_WMS(j), phi_R(j)] = Phi_VAR(MI_calc,Gs);
            red_phi(j) =  phi_R(j)-phi_WMS(j);
        end
  
    end

    Phis = [phi_WMS; phi_R; red_phi];
%     fprintf("done! number of excluded values were " + sprintf('%01d',nerr)+"\n");
end


%%%% Utility functions %%%%

% check MI from the optimisation of PID
function [flag] = check_MI_PID(PID, MI)
    flag = 1;
    MI_calc = sum(PID);
        
    if(abs(MI_calc-MI)>MI/1e4)
        disp("Mismatch between true MI and MI calculated during the optimisation!")
        disp(abs(MI_calc-MI)); %disp(MI); 
        flag = 0;
    end
end

% check MI from the optimisation of PhiID
function [flag] = check_MI_PhiID(atoms, MI)
    flag = 1;
    fn = fieldnames(atoms); MI_calc=0;
    for i = 1:length(fn)
        MI_calc = MI_calc+atoms.(fn{i});
    end

    if(abs(MI_calc-MI)>MI/1e4) 
        disp("Mismatch between true MI and MI calculated during the optimisation!")
        disp(abs(MI_calc-MI)); %disp(MI); 
        flag=0; 
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