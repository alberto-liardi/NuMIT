%% SaveQuantPID
% 
% Saves the quantiles of PID atoms to a .csv file.
% 
% Description:
%   This function takes in vectors of quantiles for PID atoms and the corresponding 
%   mutual information (MI), and saves them to a .csv file. 
%   The file name is determined by the input path and the model order 
%   (k=1 for VAR(1) models, k=2 for generic VAR(p) models).
% 
% Inputs:
%   UnX   - [L] array of quantile values for Unique X atoms.
%   UnY   - [L] array of quantile values for Unique Y atoms.
%   Red   - [L] array of quantile values for Redundancy atoms.
%   Syn   - [L] array of quantile values for Synergy atoms.
%   MI    - [L] array of total mutual information values corresponding to the PID atoms.
%   path  - String specifying the destination path for the .csv file.
%   k     - Model order flag (1 for VAR(1) model, 2 for generic VAR(p) model).
% 
% Outputs:
%   None. The function does not return any outputs. 
%   It saves a .csv file at the specified path with the name "Quantiles_PID_1.csv" if `k=1` or "Quantiles_PID_p.csv" if `k=2`.
% 
% Notes:
%   - The .csv file includes the saved variables "UnX", "UnY", "Red", "Syn", and "MI".
% 
% Alberto Liardi, 2024


function SaveQuantPID(UnX, UnY, Red, Syn, MI, path, k)

    if ~exist('k','var') || (k~=1 && k~=2), error("Insert valid parameter"); end
    
    app = ["_1","_p"];
    
    assert(length(MI)==length(UnX) && length(MI)==length(UnY) ...
       && length(MI)==length(Red) && length(MI)==length(Syn), ...
       "Dimensions of PID atoms do not match!");
    
    qtab = array2table([UnX; UnY; Red; Syn]');
    qtab.Properties.VariableNames(1:4) = {'Unique_Information_X', ...
                        'Unique_Information_Y', 'Redundancy', 'Synergy'};
    writetable(qtab, path+"/Quantiles_PID"+app(k)+".csv"); 

end