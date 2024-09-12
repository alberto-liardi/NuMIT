%% Null_model_Gauss
% 
% Computes the null distribution of PID atoms for Gaussian systems.
% 
% Description:
%   This function computes the null distribution of PID atoms for Gaussian systems of the form 
%   T = A*S+e, for a specified value of mutual information (MI). The function generates 
%   random Gaussian systems with a given MI and computes PID atoms by using different redundancy 
%   functions (MMI, DEP, or CCS). It returns matrices of the computed PID atoms across multiple 
%   runs, representing the null distribution.
% 
% Inputs:
%   MI      - Target mutual information value for the Gaussian system.
%   Sources - Number of sources in the system (default is 2).
%   Targets - Number of targets in the system (default is 1).
%   N_runs  - Number of runs for the simulation (default is 10,000).
%   red_fun - Redundancy function to use: "MMI" (default), "DEP", "CCS", or "all" (for all functions).
% 
% Outputs:
%   PIDs    - [4, N_runs, 3] matrix where each slice represents PID atoms computed using different
%             redundancy functions:
%             [:,:,1] MMI definition
%             [:,:,2] DEP definition
%             [:,:,3] CCS definition
%             Each of the N_runs column contains PID atoms in the order:
%             [Unique X; Unique Y; Redundancy; Synergy]
% 
% Notes:
%   - The function verifies the results and skips runs if the MI value is not reached or if errors
%     occur during calculations.
%   - The output matrix has three slices, even if only one redundancy function is selected.
% 
% Alberto Liardi, 2024

function PIDs = Null_model_Gauss(MI,Sources,Targets,N_runs,red_fun)

    if (~exist('Sources','var') || isempty(Sources)), Sources = 2; end
    if ~exist('Targets','var') || isempty(Targets), Targets = 1; end
    if ~exist('N_runs','var') || isempty(N_runs), N_runs = 1e4; end
    if ~exist('red_fun','var') || isempty(red_fun), red_fun = "MMI"; end

    Syn = NaN(3,N_runs); Red = NaN(3,N_runs);
    UnX = NaN(3,N_runs); UnY = NaN(3,N_runs);
    PIDs = NaN(4,N_runs,3);
    
    nerr = 0;
    S = Sources;
    T = Targets;
        
    for j = 1:N_runs
        
        % if(mod(j,N_runs/10)==0), disp(j); end
        
        % generate random Gaussian system with specific MI
        [Sigma_full, check] = Gaussian_from_MI(MI,S,T);
    
        if(check==1)
            disp("problem with MI optimisation"); 
            nerr = nerr+1;
            continue
        end

        
        % now calculate PID atoms using {red_fun} definition
        if red_fun=="MMI" || red_fun=="all"
            % MMI
            [UnX(1,j), UnY(1,j), Red(1,j), Syn(1,j)] = PID_MMI_Gaussian(Sigma_full,S,T);
        end
        if red_fun=="DEP" || red_fun=="all"
            % DEP
            [UnX(2,j), UnY(2,j), Red(2,j), Syn(2,j)] = PID_DEP_Gaussian(Sigma_full,S,T);
        end
        if red_fun=="CCS" || red_fun=="all"
            % CCS 
            [UnX(3,j), UnY(3,j), Red(3,j), Syn(3,j)] = PID_CCS_Gaussian(Sigma_full,S,T);
        end
    
        % check consistency of PID atoms
        try
            PIDcheck([UnX(1,j), UnY(1,j), Red(1,j), Syn(1,j)]);
        catch
            nerr = nerr+1;
            continue
        end
        assert(isreal(Syn(1,j)), "not real");
    
    end

    for n = 1:3
        PIDs(:,:,n) = [UnX(n,:); UnY(n,:); Red(n,:); Syn(n,:)];
    end
end
