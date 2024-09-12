%% PhiID_VAR_calculator
%
% Computes the PhiID atoms for a VAR(p) model with multivariate sources.
% 
% Description:
%   This function calculates the PhiID (Partial Information Decomposition) atoms for a VAR(p) model
%   with multivariate sources. The PhiID atoms represent various information-theoretic quantities
%   such as redundancy, synergy, and unique information within the system.
%   It returns the PhiID atoms in a struct and the autocovariance matrices of the VAR model.
% 
% Inputs:
%   p       - Order of the VAR(p) model.
%   A       - [p] cell array of [L1+L2, L1+L2] matrix coefficients of the VAR model.
%   V       - [L1+L2, L1+L2] residual covariance matrix.
%   L1      - Length of the first source.
%   L2      - Length of the second source.
% 
% Outputs:
%   atoms   - Struct containing PhiID atoms with the following fields:
%             rtr - {1}{2}->{1}{2} atom
%             rtx - {1}{2}->{1} atom
%             rty - {1}{2}->{2} atom
%             rts - {1}{2}->{12} atom
%             xtr - {1}->{1}{2} atom
%             xtx - {1}->{1} atom
%             xty - {1}->{2} atom
%             xts - {1}->{12} atom
%             ytr - {2}->{1}{2} atom
%             ytx - {2}->{1} atom
%             yty - {2}->{2} atom
%             yts - {2}->{12} atom
%             str - {12}->{1}{2} atom
%             stx - {12}->{1} atom
%             sty - {12}->{2} atom
%             sts - {12}->{12} atom
%   Gamma   - {p+1} cell array of autocovariance matrices of the VAR model.
% 
% Notes:
%   - The function currently supports only VAR(1) models (p=1).
%   - Ensure that the coefficient matrices are provided as a cell array.
%   - The function solves the Lyapunov equation to obtain the autocovariance matrices.
%   - Mutual information values and PhiID atoms are computed in bits (log base 2).
% 
% Alberto Liardi, 2024


function [atoms, Gamma] = PhiID_VAR_calculator(p,A,V,L1,L2)

    L = L1 + L2; % length of the target
    
    if ~exist('V','var'), V = eye(L); end
    assert(iscell(A), "Coefficient matrices must be in a 1xp cell");
    % currently working only for p=1 (VAR(1) systems)
    assert(p==1, "Enter a valid model order"); 
    assert(length(A)==p, "Model order and number of evolution matrices are not the same");
    assert(size(A{1},1) == L, "Sources and evolution matrix dimensions are not compatible");
    
    
    % transform the cell array A into a 3D matrix
    A = cat(3, A{:});
    
    % solving the Lyapunov equation
    G = var_to_autocov(A,V,p);
    
    % rewrite it with a cell array
    Gamma = cell(1,p+1);
    for s = 1:p+1
        Gamma{s} = G(:,:,s);
    %     disp(Gamma{s});
    end
    
    
    % let's now construct the covariance matrices
    
    % T = (a, b) are the targets
    % S = (x, y) are the sources
    
    % full covariance matrix and its entropy
    Sigma_full = cell2mat({Gamma{1}, Gamma{2}; Gamma{2}', Gamma{1}});
    H_xyab =h(Sigma_full);
    
    % source covariance and its entropy
    Sigma_xy = Gamma{1}; H_xy = h(Sigma_xy);
    
    % target covariance and its entropy
    Sigma_ab = Gamma{1}; H_ab = h(Sigma_ab);
    
    % single sources and targets (x-a and y-b are equal due to stationarity)
    Sigma_x = Gamma{1}(1:L1,1:L1); H_x = h(Sigma_x);
    Sigma_a = Sigma_x; H_a = h(Sigma_a);
    Sigma_y = Gamma{1}(1+L1:L1+L2,1+L1:L1+L2); H_y = h(Sigma_y);
    Sigma_b = Sigma_y; H_b = h(Sigma_b);
    
    % combinations of sources and targets
    Sigma_xab = cell2mat({Gamma{1}(1:L1,1:L1), Gamma{2}(:,1:L1)'; ...
                          Gamma{2}(:,1:L1), Gamma{1}}); 
    H_xab = h(Sigma_xab);
    
    Sigma_yab = cell2mat({Gamma{1}(1+L1:L1+L2,1+L1:L1+L2),Gamma{2}(:,1+L1:L1+L2)'; ...
                          Gamma{2}(:,1+L1:L1+L2),Gamma{1}});
    H_yab = h(Sigma_yab);
     
    Sigma_xya = cell2mat({Gamma{1},Gamma{2}(1:L1,:)'; Gamma{2}(1:L1,:),Gamma{1}(1:L1,1:L1)});
    H_xya = h(Sigma_xya);
    
    Sigma_xyb = cell2mat({Gamma{1},Gamma{2}(1+L1:L1+L2,:)'; ...
                         Gamma{2}(1+L1:L1+L2,:),Gamma{1}(1+L1:L1+L2,1+L1:L1+L2)});
    H_xyb = h(Sigma_xyb);
    
    % mixed pairs of sources and targets 
    Sigma_xa = cell2mat({Gamma{1}(1:L1,1:L1),Gamma{2}(1:L1,1:L1)'; ...
                         Gamma{2}(1:L1,1:L1),Gamma{1}(1:L1,1:L1)});
    H_xa = h(Sigma_xa);
    
    Sigma_xb = cell2mat({Gamma{1}(1:L1,1:L1),Gamma{2}(1+L1:L1+L2,1:L1)'; ...
                         Gamma{2}(1+L1:L1+L2,1:L1),Gamma{1}(1+L1:L1+L2,1+L1:L1+L2)});
    H_xb = h(Sigma_xb);
    
    Sigma_ya = cell2mat({Gamma{1}(1+L1:L1+L2,1+L1:L1+L2),Gamma{2}(1:L1,1+L1:L1+L2)'; ...
                         Gamma{2}(1:L1,1+L1:L1+L2),Gamma{1}(1:L1,1:L1)});
    H_ya = h(Sigma_ya);
     
    Sigma_yb = cell2mat({Gamma{1}(1+L1:L1+L2,1+L1:L1+L2),Gamma{2}(1+L1:L1+L2,1+L1:L1+L2)'; ...
                         Gamma{2}(1+L1:L1+L2,1+L1:L1+L2),Gamma{1}(1+L1:L1+L2,1+L1:L1+L2)});
    H_yb = h(Sigma_yb);
    
    % checks
    assert(all(all(Sigma_x-Sigma_x'<1e-8))); 
    assert(all(all(Sigma_y - Sigma_y' < 1e-8)));
    assert(all(all(Sigma_a - Sigma_a' < 1e-8)));
    assert(all(all(Sigma_b - Sigma_b' < 1e-8)));
    assert(all(all(Sigma_xa - Sigma_xa' < 1e-8)));
    assert(all(all(Sigma_xb - Sigma_xb' < 1e-8)));
    assert(all(all(Sigma_ya - Sigma_ya' < 1e-8)));
    assert(all(all(Sigma_yb - Sigma_yb' < 1e-8)));
    assert(all(all(Sigma_xab - Sigma_xab' < 1e-8)));
    assert(all(all(Sigma_yab - Sigma_yab' < 1e-8)));
    assert(all(all(Sigma_xya - Sigma_xya' < 1e-8)));
    assert(all(all(Sigma_xyb - Sigma_xyb' < 1e-8)));
    assert(all(all(Sigma_full - Sigma_full' < 1e-8)));
    assert(all(all(Sigma_xy - Sigma_xy' < 1e-8)));
    assert(all(all(Sigma_ab - Sigma_ab' < 1e-8)));
     
    chol(Sigma_x); chol(Sigma_y); chol(Sigma_a); chol(Sigma_b);
    chol(Sigma_xa); chol(Sigma_xb); chol(Sigma_ya); chol(Sigma_yb);
    chol(Sigma_xab); chol(Sigma_yab); chol(Sigma_xya); chol(Sigma_xyb);
    chol(Sigma_full); chol(Sigma_xy); chol(Sigma_ab);
    
    % now calculate the mutual information
    MI_xytab = (H_xy + H_ab- H_xyab) / log(2);
    
    MI_xta = (H_x + H_a - H_xa) / log(2);
    MI_yta = (H_y + H_a - H_ya) / log(2);
    MI_xtb = (H_x + H_b - H_xb) / log(2);
    MI_ytb = (H_y + H_b - H_yb) / log(2);
    
    MI_xyta = (H_xy + H_a - H_xya) / log(2);
    MI_xytb = (H_xy + H_b - H_xyb) / log(2);
    MI_xtab = (H_x + H_ab - H_xab) / log(2);
    MI_ytab = (H_y + H_ab - H_yab) / log(2);
    
    Rxyta = Red(MI_xta, MI_yta);
    Rxytb = Red(MI_xtb, MI_ytb);
    Rabtx = Red(MI_xta, MI_xtb);
    Rabty = Red(MI_yta, MI_ytb);
    Rxytab = Red(MI_xtab, MI_ytab);
    Rabtxy = Red(MI_xyta, MI_xytb);
    
    doubleR = doubleRed(MI_xta, MI_xtb, MI_yta, MI_ytb);
    
    
    % solve the linear system of equations
    MIs = [doubleR Rxyta Rxytb Rxytab Rabtx Rabty Rabtxy ...
            MI_xta MI_xtb MI_yta MI_ytb MI_xyta MI_xytb MI_xtab MI_ytab MI_xytab];
    
    Mat = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; % doubleR
         1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0; % Rxyta
         1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0; % Rxytb
         1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0; % Rxytab
         1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0; % Rabtx
         1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0; % Rabty
         1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0; % Rabtxy
         1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0; % MI_xta
         1 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0; % MI_xtb
         1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0; % MI_yta
         1 0 1 0 0 0 0 0 1 0 1 0 0 0 0 0; % MI_ytb
         1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0; % MI_xyta
         1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0; % MI_xytb
         1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0; % MI_xtab
         1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0; % MI_ytab
         1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; % MI_xytab
    
    lattice = linsolve(Mat, MIs');
    
    % Sort the results and return
    atoms = [];
    atoms.rtr = lattice(1,:);
    atoms.rtx = lattice(2,:);
    atoms.rty = lattice(3,:);
    atoms.rts = lattice(4,:);
    atoms.xtr = lattice(5,:);
    atoms.xtx = lattice(6,:);
    atoms.xty = lattice(7,:);
    atoms.xts = lattice(8,:);
    atoms.ytr = lattice(9,:);
    atoms.ytx = lattice(10,:);
    atoms.yty = lattice(11,:);
    atoms.yts = lattice(12,:);
    atoms.str = lattice(13,:);
    atoms.stx = lattice(14,:);
    atoms.sty = lattice(15,:);
    atoms.sts = lattice(16,:);
    
end

function entropy = h(S)
    entropy = 0.5*logdet(S);
end

function I = Red(I1, I2)
    I = min([I1, I2]);
end

function I = doubleRed(I1, I2, I3, I4)
    I = min([I1, I2, I3, I4]);
end