%% PID_VAR_calculator
%
% Computes the PID atoms for a VAR(p) model with multivariate sources.
% 
% Description:
%   This function calculates the Partial Information Decomposition (PID) atoms for a VAR(p) model
%   with multivariate sources. The function supports different redundancy functions (MMI, DEP, CCS)
%   and allows specifying the target of PID calculations. It returns a matrix containing the PID atoms 
%   for each redundancy definition and the autocovariance matrices of the VAR model.
% 
% Inputs:
%   p          - Order of the VAR(p) model.
%   A          - {p} cell array of [L1+L2, L1+L2] matrix coefficients of the VAR model.
%   V          - [L1+L2, L1+L2] residual covariance matrix.
%   L1         - Length of the first source.
%   L2         - Length of the second source.
%   red_fun    - Redundancy functions to use: "MMI" (default), "DEP", "CCS", or "all" (for all functions).
% 
% Outputs:
%   PID        - [3,4] matrix where the first row contains MMI-PID atoms, the second row contains
%             DEP-PID atoms, and the third row contains CCS-PID atoms. The atoms are ordered as:
%             [UniqueX, UniqueY, Redundancy, Synergy]. If a redundancy function is not used, zeros 
%             are returned for its row.
%   Gamma      - {p+1} cell array of autocovariance matrices of the VAR model.
% 
% Alberto Liardi, 2024


function [PID, Gamma] = PID_VAR_calculator(p,A,V,L1,L2,red_fun)

L = L1 + L2; % length of the target

if ~exist('V','var'), V = eye(L); end
assert(iscell(A), "Coefficient matrices must be in a 1xp cell");
assert(0<p, "Enter a valid model order"); 
assert(length(A)==p, "Model order and number of evolution matrices are not the same");
assert(size(A{1},1) == L, "Sources and evolution matrix dimensions are not compatible");
if ~exist('red_fun','var') || isempty(red_fun), red_fun = ["MMI"]; end

% creating the matrices for the Lyapunov equation

B = zeros(L,L,p);
for s = 1:p
    B(:,:,s) = A{s};
end
G = var_to_autocov(B,V,p);
Gamma = cell(1,p+1);
for s = 1:p+1
    Gamma{s} = G(:,:,s);
end

% T = (T1, T2) is the target
% S = (S1, S2) is the source
% TS all of them (target(s) + sources)

sigma_cov_TS = cell(p+1,p+1);
sigma_cov_TS(:) = {zeros(L,L)};
for s = 1:1:p+1
    for t = s:1:p+1
        sigma_cov_TS{s,t} = Gamma{t-s+1};
    end
end
% going from the upper diagonal to the full matrix 
sigma_cov_TS = cell2mat(sigma_cov_TS);
sigma_cov_TS = sigma_cov_TS - tril(sigma_cov_TS,-1);
sigma_cov_tot_TS= triu(sigma_cov_TS,1) + sigma_cov_TS';

sigma_cov_tot_S = sigma_cov_tot_TS(1:end-L, 1:end-L);
%     disp(sigma_cov_tot_S);

sigma_cov_tot_T = Gamma{1};
%     disp(sigma_cov_tot_T);

H_T = 0.5*logdet(sigma_cov_tot_T);
H_S = 0.5*logdet(sigma_cov_tot_S);
H_TS = 0.5*logdet(sigma_cov_tot_TS);

% mutual information
MI = (H_T + H_S - H_TS) / log(2);    

% some more covariance matrices for the marginals
sigma_T1S1 = zeros(L1*(p+1),L1*(p+1));
sigma_T2S2 = zeros(L2*(p+1),L2*(p+1));
for s = 1:1:p+1
    for t = s:1:p+1
        sigma_T1S1(s*L1-L1+1:s*L1,t*L1-L1+1:t*L1) = Gamma{t-s+1}(1:L1,1:L1);
        sigma_T2S2(s*L2-L2+1:s*L2,t*L2-L2+1:t*L2) = Gamma{t-s+1}(L1+1:end,L1+1:end);
    end
end

% concatenating matrices to obtain the full covariance matrix (_tot)
sigma_T1S1_tot = sigma_T1S1' + triu(sigma_T1S1,1) - tril(sigma_T1S1,-1)';
sigma_T2S2_tot = sigma_T2S2' + triu(sigma_T2S2,1) - tril(sigma_T2S2,-1)';

row_T2 = zeros(L2,L1*(p+1));
row_T1 = zeros(L1,L2*(p+1));
for l = 1:1:p+1
    row_T2(1:L2,l*L1-L1+1:l*L1) = Gamma{l}(L1+1:end,1:L1);
    row_T1(1:L1,l*L2-L2+1:l*L2) = Gamma{l}(1:L1,L1+1:end);
end

addT2 = Gamma{1}(L1+1:end,L1+1:end);
addT1 = Gamma{1}(1:L1,1:L1);
for l = 1:1:p+1
    row_T2(1:L2,l*L1-L1+1:l*L1) = Gamma{l}(L1+1:L1+L2,1:L1);
    row_T1(1:L1,l*L2-L2+1:l*L2) = Gamma{l}(1:L1,L1+1:L1+L2);
end

sigma_TS1 = vertcat( row_T2, sigma_T1S1_tot );
column = vertcat( addT2, row_T2' );
sigma_TS1 = horzcat( column, sigma_TS1);

sigma_TS2 = vertcat( row_T1, sigma_T2S2_tot );
column = vertcat( addT1, row_T1' );
sigma_TS2 = horzcat( column, sigma_TS2);

sigma_S1 = sigma_T1S1_tot(L1+1:end,L1+1:end);
sigma_S2 = sigma_T2S2_tot(L2+1:end,L2+1:end);

% calculating I(T,S1):
H_T = 0.5*logdet(sigma_cov_tot_T);
H_S = 0.5*logdet(sigma_S1);
H_TS = 0.5*logdet(sigma_TS1);

MI_TX = (H_T + H_S - H_TS) / log(2);

% calculating I(T,S2):
H_T = 0.5*logdet(sigma_cov_tot_T);
H_S = 0.5*logdet(sigma_S2);
H_TS = 0.5*logdet(sigma_TS2);

MI_TY = (H_T + H_S - H_TS) / log(2);

% need to reorder the above covariance for CCS and DEP pipelines
Covariance = sigma_cov_tot_TS;
if any(ismember(["DEP", "CCS", "all"], red_fun))
    Cov_reord = CovReord(Covariance,L1,L2,p); 
end

% now calculating PID quantities using the 'red_fun' definitions
PID = zeros(3,4);
if any(ismember(["MMI", "all"], red_fun))
    % MMI
    R = min(MI_TX, MI_TY); 
    U_X = MI_TX - R;
    U_Y = MI_TY - R;
    S = MI - U_Y - U_X - R;
    
    print = "MMI Redundancy is: " + R + ...
        "\nUnique information of X is: " + U_X + ...
        "\nUnique information of Y is: " + U_Y + ...
        "\nSynergy is: " + S + "\n\n";
%         fprintf(print); 

    PID(1,:) = [U_X, U_Y, R, S];
end        
if any(ismember(["DEP", "all"], red_fun))
    % DEP
    PI = num2cell(calc_pi_Idep_mvn(Cov_reord, [p*L1 p*L2 L]));
    [R, U_X, U_Y, S] = PI{:};

    print = "DEP Redundancy is: " + R + ...
        "\nUnique information of X is: " + U_X + ...
        "\nUnique information of Y is: " + U_Y + ...
        "\nSynergy is: " + S + "\n\n";
%         fprintf(print); 

    PID(2,:) = [U_X, U_Y, R, S];
end
if any(ismember(["CCS", "all"], red_fun))
    % CCS 
    try 
        PI = num2cell(calc_pi_mvn(lattice2d(), Cov_reord, ...
                                [p*L1 p*L2 L], @Iccs_mvn_P2).PI);
        [R, U_X, U_Y, S] = PI{:};
    
        print = "CCS Redundancy is: " + R + ...
            "\nUnique information of X is: " + U_X + ...
            "\nUnique information of Y is: " + U_Y + ...
            "\nSynergy is: " + S + "\n\n";
        % fprintf(print); 

        PID(3,:) = [U_X, U_Y, R, S];

    catch 
        disp("CCS calculation failed, skipping...");
        PID(3,:) = 0;
    end
end

end