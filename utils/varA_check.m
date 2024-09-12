%% varA_check
% 
% Checks if the coefficient matrices of a VAR model have too small spectral radius.
% 
% Description:
%   This function checks whether the coefficient matrices of a VAR model have spectral radius that are 
%   below the threshold value (1e-18). 
%   It returns a boolean indicating whether the coefficient matrices pass the check.
% 
% Inputs:
%   A    - Cell array of coefficient matrices from the VAR model.
% 
% Outputs:
%   bool - Returns 1 if no values smaller than 1e-18 are found, and 0 otherwise.
% 
% 
% Alberto Liardi, 2024


function bool = varA_check(A)

    if isempty(find(abs(A)<1e-18,1))
        bool = 1;
    else 
        bool = 0;
    end

end



    