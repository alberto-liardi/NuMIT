%% LoadPhiID
% 
% Loads PhiID atoms from a .csv file and computes the total mutual information (MI).
% 
% Description:
%   This function reads PhiID atoms from a specified file and converts the data into a cell array 
%   of structs, where each struct contains the PhiID atoms for a given sample. It also calculates 
%   the total mutual information (MI) as the sum of all PhiID atoms for each sample.
% 
% Inputs:
%   file    - Path to the file containing PhiID atom data in table format.
% 
% Outputs:
%   atoms   - Cell array of structs, where each struct contains PhiID atoms for one sample.
%   MI      - [1xL] array where each element is the total mutual information for the corresponding sample.
% 
% Alberto Liardi, 2024


function [atoms, MI] = LoadPhiID(file)

    if ~isfile(file), error(file+" directory has not been found!"); end

    T = readtable(file);
    Tstr= table2struct(T);
    
    L = length(Tstr);
    % convert array of struct into cell of structs
    atoms = cell(1, L);
    for j = 1:length(Tstr)
        atoms{j} = Tstr(j);
    end
    
    % calculate total TDMI (sum of all atoms)
    MI = zeros(1, L);
    fn = fieldnames(atoms{1});
    for j = 1:L  
        for i = 1:length(fn)
            MI(j) = MI(j)+atoms{j}.(fn{i});
        end
    end

end