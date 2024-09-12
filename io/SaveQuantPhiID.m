%% SaveQuantPhiID
% 
% Saves PhiID atoms to a .csv file.
% 
% Description:
%   This function converts a matrix of quantile PhiID atoms into a table format and saves it as a .csv file. 
% 
% Inputs:
%   qatoms  - [16, L] matrix of quantile PhiID atoms, where each row corresponds to a specific PhiID atom type:
%             * rtr - rtr atom
%             * rtx - rtx atom
%             * rty - rty atom
%             * rts - rts atom
%             * xtr - xtr atom
%             * xtx - xtx atom
%             * xty - xty atom
%             * xts - xts atom
%             * ytr - ytr atom
%             * ytx - ytx atom
%             * yty - yty atom
%             * yts - yts atom
%             * str - str atom
%             * stx - stx atom
%             * sty - sty atom
%             * sts - sts atom
%   outpath - Path to the directory where the .csv file will be saved.
% 
% Outputs:
%   None. The function does not return any outputs. 
%   It saves a .csv file named `Quantiles_PhiID.csv` in the specified directory (`outpath`). 
%   The file contains the quantile PhiID atoms formatted in a table.
% 
% Alberto Liardi, 2024


function [] = SaveQuantPhiID(qatoms, outpath)

    % convert matrix in a table and save it
    qtab = array2table(qatoms');
    qtab.Properties.VariableNames(1:16) = {'rtr', 'rtx', 'rty', 'rts', 'xtr', ...
                                           'xtx', 'xty', 'xts', 'ytr', 'ytx', ...
                                           'yty', 'yts', 'str', 'stx', 'sty', 'sts'};
    writetable(qtab, outpath+'/Quantiles_PhiID.csv'); 

end