%% Null_model_Discrete
% 
% Generates the null distribution of PID atoms for discrete systems.
% 
% Description:
%   This function computes the null distribution of PID atoms for discrete systems of the form T=f(X,Y)+pe,
%   for a specified value of mutual information (MI). The function generates random discrete systems 
%   with different 2-bit logic gates and porbabilities from a Dirichlet distribution, then calculates 
%   PID atoms by optimizing the noise parameter pe to match the specified MI.
%   It returns matrices of the computed PID atoms across multiple runs, i.e. the null distribution.
% 
% Inputs:
%   MI      - Target mutual information value for the discrete system.
%   N_runs  - Number of runs for the simulation (default is 70,000).
% 
% Outputs:
%   PIDs    - [4, N_runs] matrix where each of the N_runs columns contains PID atoms in the order:
%             [Unique X; Unique Y; Redundancy; Synergy].
% 
% Notes:
%   - The function verifies the results and skips steps if the MI value is not reached or if the optimization fails.
% 
% Alberto Liardi, 2024

function PIDs = Null_model_Discrete(MI,N_runs)

    if ~exist('N_runs','var') || isempty(N_runs), N_runs = 7e4; end

    nerr = 0;
    
    % discrete systems have a maximum of MI at -log(0.5):
    assert(MI<-log(0.5), ...
           "Value of Mutual Information is incompatible with the system!")
    
    Syn = NaN(1,N_runs); Red = NaN(1,N_runs);
    UnX = NaN(1,N_runs); UnY = NaN(1,N_runs);
    
    % sample the parameters for the dirichlet from a uniform distribution
    a = repmat(unifrnd(0.1, 3),1,4);

    % set of indices for the loop
    runarr = 1:floor(N_runs/7);
    
    for gate = ["XOR","XOR2","XOR3","OR","OR2","OR3","OR4"]

        fprintf("doing gate %s!\n", gate);

        if gate == "XOR"
            boolfun = @XORgate;
            js = runarr;
        elseif gate == "XOR2"
            boolfun = @XOR2gate;
            js = floor(N_runs/7)+[runarr];
        elseif gate == "XOR3"
            boolfun = @XOR3gate;
            js = floor(2*N_runs/7)+[runarr];
        elseif gate == "OR"
            boolfun = @ORgate;
            js = floor(3*N_runs/7)+[runarr];
        elseif gate == "OR2"
            boolfun = @OR2gate;
            js = floor(4*N_runs/7)+[runarr];
        elseif gate == "OR3"
            boolfun = @OR3gate;
            js = floor(5*N_runs/7)+[runarr];
        elseif gate == "OR4"
            boolfun = @OR4gate;
            js = floor(6*N_runs/7)+[runarr];
        end
    
        for j = js

%             if mod(j,N_runs/10)==0, disp(j), end

            % sample the probabilities of the Sources S from the dirichlet
            ps = drchrnd(1,a);

            % generate random Discrete system with specific MI
            [pe, check] = Discrete_from_MI(MI,ps,boolfun);
        
            if(check==1)
%                 disp("problem with MI optimisation"); 
                nerr = nerr+1;
                continue
            end

            [UnX(j), UnY(j), Red(j), Syn(j), MI_calc] = PID_MMI_Discrete(ps, pe, gate);
            
            if(abs(MI_calc-MI)>MI/1e4), nerr=nerr+1; 
                disp("Mismatch between true MI and MI calculated during the optimisation!")
                disp(abs(MI_calc-MI)); 
                continue; 
            end
        
            PIDcheck([UnX(j), UnY(j), Red(j), Syn(j)]);
    
        end
    end

    PIDs = [UnX; UnY; Red; Syn];

end

%%% gate functions %%%
function [pz0,pz1] = XORgate(ps)
    pz0 = ps(1)+ps(4);
    pz1 = 1-pz0;
end

function [pz0,pz1] = XOR2gate(ps)
    pz0 = ps(1)+ps(2);
    pz1 = 1-pz0;
end

function [pz0,pz1] = XOR3gate(ps)
    pz0 = ps(1)+ps(3);
    pz1 = 1-pz0;
end

function [pz0,pz1] = ORgate(ps)
    pz0 = ps(1);
    pz1 = 1-ps(1);
end

function [pz0,pz1] = OR2gate(ps)
    pz0 = ps(2);
    pz1 = 1-pz0;
end

function [pz0,pz1] = OR3gate(ps)
    pz0 = ps(3);
    pz1 = 1-pz0;
end

function [pz0,pz1] = OR4gate(ps)
    pz0 = ps(4);
    pz1 = 1-pz0;
end


function r = drchrnd(n,a)
    p = length(a);
    r = gamrnd(repmat(a,n,1),1,n,p);
    r = r ./ repmat(sum(r,2),1,p);
end