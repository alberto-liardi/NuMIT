%% cellstruct2mat
% 
% Converts a cell array of structs into a 2D matrix.
% 
% Description:
%   This function takes a cell array of structs and converts it into a 2D matrix. Each struct in the cell array 
%   is converted into a column of the matrix, where the number of rows corresponds to the number of fields in 
%   each struct. The resulting matrix has dimensions [number of fields x length of the cell array].
% 
% Inputs:
%   s   - Cell array of structs, where each struct contains the same fields.
% 
% Outputs:
%   mat - [M x N] matrix, where M is the number of fields in the struct and N is the length of the cell array.
%         Each column of the matrix corresponds to the values of a struct in the cell array.
% 
% Alberto Liardi, 2024


function mat = cellstruct2mat(s)

    mat = NaN(length(fieldnames(s{1})), length(s));
    for i = 1:length(s)
        % each cell becomes an array (of length = length(struct))
        mat(:,i) = cell2mat(struct2cell(s{i}));
    end
end