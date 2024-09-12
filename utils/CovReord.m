%% CovReord
% 
% Reorder a covariance matrix from a specific block structure to a new structure.
% 
% Description:
%   This function reorders a covariance matrix originally structured as follows:
%   [T (target), S1_1 (first source first lag), S2_1 (second source first lag), ...,
%   ..., S1_p (first source p-th lag), S2_p (second source p-th lag)],
%   to a matrix where the sources are grouped together and the target is at the end:
%   [S1 (first source all lags), S2 (second source all lags), T (target)].
% 
% Inputs:
%   C     - [n,n] covariance matrix with the original block structure.
%   L1    - Number of variables in the first source.
%   L2    - Number of variables in the second source.
%   p     - Number of lags for each source.
% 
% Outputs:
%   C     - [n,n] reordered covariance matrix with the new block structure: S1, S2, T.
% 
% Alberto Liardi, 2024


function C = CovReord(C,L1,L2,p)
    
    % control the input dimensions
    assert( size(C,1) == size(C,2), "Covariance matrix is not a square matrix!");
    assert( size(C,1) == (L1+L2)*(p+1), ...
            "Dimension of covariance matrix and sources-targets are not compatible");
    
    vs = [L1*p L2*p L1+L2];
    L = L1+L2;
    
    % indices
    odd = 1:2:sum(vs);
    even = 2:2:sum(vs);
    xx = [1:L1*(p+1)];
    yy = [L1*(p+1)+1:L1*(p+1)+L2*(p+1)];
    
    
    % first I need to move S1 to the left, and S2 to the right
    % target T is still at the beginning, will move it later
    
    new_odd = zeros(L1,length(odd));
    new_even = zeros(L2,length(even));
    
    for j = 1:length(odd)
        new_odd(:,j) = (odd(j)-1)/2*L + [1:L1];
    end
    for j = 1:length(even)    
        new_even(:,j) = (even(j)-2)/2*L + L1 + [1:L2];
    end
    
    odd = reshape(new_odd, 1,size(new_odd,1)*size(new_odd,2));
    even = reshape(new_even, 1,size(new_even,1)*size(new_even,2));
    
    odd = odd(1:L1*(p+1));
    even = even(1:L2*(p+1));
    
    Cov = C; Cov2 = C; 
    
    % moving rows
    Cov2(xx,:) = C(odd, :);
    Cov2(yy,:) = C(even, :);
    % moving columns
    Cov(:, xx) = Cov2(:, odd);
    Cov(:, yy) = Cov2(:, even);
    
    % disp("T1 - S1 - T2 - S2");
    

    % now need to reorder sources among themselves (only if multivariate)
    % (have a format like S1_X S1_Y ... S1_Z, i.e. must be separated)
    
    % Source 1
    Cov3 = Cov;
    l = 0:L1;
    for j = 1:L1          
        pos = l(j)*(p+1)+1:1:(l(j)+1)*(p+1);
        index = [j:L1:L1*(p+1)];
        % moving rows
        Cov3(pos,:) = Cov(index,:);
    end
    
    % disp("move S1 rows");
    
    for j = 1:L1   
        pos = l(j)*(p+1)+1:1:(l(j)+1)*(p+1);
        index = [j:L1:L1*(p+1)];
        % moving columns
        Cov(:,pos) = Cov3(:, index);
    end
    
    Cov(:,L1*(p+1)+1:end) = Cov3(:, L1*(p+1)+1:end);
    
    % disp("move S1 col");
    
    % now need to do the same for the second source S2

    % Source 2
    Cov2 = Cov;
    l = 0:L2;
    for j = 1:L2           
        pos = L1*(p+1)+l(j)*(p+1)+1:1:(l(j)+1)*(p+1)+L1*(p+1);
        index = [L1*(p+1)+j:L2:L*(p+1)];
        % moving rows
        Cov2(pos,:) = Cov(index,:);
    end
    
    % disp("move S2 rows"); Cov2
    
    for j = 1:L2   
        pos = L1*(p+1)+l(j)*(p+1)+1:1:(l(j)+1)*(p+1)+L1*(p+1);
        index = [L1*(p+1)+j:L2:L*(p+1)];
        % moving columns
        Cov(:,pos) = Cov2(:, index);
    end
    
    Cov(:,1:L1*(p+1)) = Cov2(:,1:L1*(p+1));
    
    % disp("move S2 cols"); Cov
    
    C = Cov;
    
    % moving target at the end and reordering (target is the joint future)
    
    T1 = 1:p+1:L1*(p+1);
    T2 = L1*(p+1)+1:p+1:L*(p+1);
    
    row_T1 = Cov(T1,:); % take T1 rows
    row_T2 = Cov(T2,:); % take T2 rows
    
    dim = L*(p+1);
    
    % all targets
    T = horzcat(T1, T2);
    
    % avoiding targets
    nT = setdiff([1:dim], T);
    
    row_T1 = horzcat( row_T1(:,nT), row_T1(:,T));
    row_T2 = horzcat( row_T2(:,nT), row_T2(:,T));
    
    C = C(nT, nT);
    
    C = horzcat(C, row_T1(:, 1:end-L)', row_T2(:, 1:end-L)');
    C = vertcat(C, row_T1, row_T2);

end