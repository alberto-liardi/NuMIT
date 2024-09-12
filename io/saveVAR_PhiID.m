%% saveVAR_PhiID
% 
% Saves PhiID atoms and related information from a VAR-PhiID procedure to files.
% 
% Description:
%   This function saves the PhiID atoms calculated from a VAR PhiID procedure to a CSV file. 
%   It also saves key parameters and subject information to separate text files for reference. 
%   The file names and paths are constructed based on subject and condition details from the input.
% 
% Inputs:
%   Information - Struct containing the following fields:
%                 * subject_info - Struct with fields:
%                     - name      : Subject name or identifier.
%                     - drug      : Name of the drug condition (for psychedelic analysis).
%                     - condition : Condition type (e.g., "task", "resting state").
%                     - ID        : Unique subject identifier.
%                 * parameters - Struct with fields:
%                     - channels  : Number of channels used in the analysis.
%                     - epochs    : Number of epochs in the data.
%                     - runs      : Number of runs in the VAR PhiID procedure.
%                 * PhiIDs    - Cell array of structs representing PhiID atoms, where each struct contains the PhiID values.
%   out         - Output directory path where files will be saved.
% 
% Outputs:
%   The function does not return any outputs. It saves the following files to the specified output directory:
%   * data.csv    - CSV file containing the PhiID atoms, with each row corresponding to a set of atoms.
%   * summary.txt - Text file summarizing the analysis parameters (channels, epochs, and number of runs).
%   * ID.txt      - Text file containing the subject's ID, name, drug, and condition information.
% 
% Alberto Liardi, 2024


function [] = saveVAR_PhiID(Information, out)

    % only useful for psychedelic analysis
    if ~isfield(Information.subject_info, 'drug'), Information.subject_info.drug=""; end
    if ~isfield(Information.subject_info, 'ID'), Information.subject_info.ID=""; end

    % creating the name of the path
    if ~isstring(Information.subject_info.name) && ~ischar(Information.subject_info.name)
        name = sprintf('%01d',Information.subject_info.name);
    else, name = Information.subject_info.name;
    end
    state = Information.subject_info.condition;
    subdir = strcat(Information.subject_info.drug,'/');
    
    path = strcat(out, subdir, name, '_c', ...
                  sprintf('%01d',Information.parameters.channels), '/', state, '/');
    if not(isfolder(path)), mkdir(path); end
    
    
    % convert the Phi_ID lattices (struct) into tables, then concatenate them
    T = struct2table(Information.PhiIDs{1});
    for i = 2:length(Information.PhiIDs)
        row = struct2table(Information.PhiIDs{i});
        T = [T;row];
    end
    
    writetable(T, path+'data.csv'); % ,'Delimiter','\t'
    
    % save some parameters
    fileID = fopen(path + 'summary.txt','w');
    fprintf(fileID,'%s \t %s \t %s \t %s \t %s \t %s \n','channels', 'epochs', 'N. runs');
    sum = [Information.parameters.channels, Information.parameters.epochs, Information.parameters.runs];
    fprintf(fileID,'%d \t %d \t %d \n', sum);
    fclose(fileID);
           
    fileID = fopen(path + 'ID.txt','w');
    fprintf(fileID,'%s \t %s \t %s \t %s \t \n','ID', 'Subject', 'Drug', 'Condition');
    sum = [Information.subject_info.ID, name, ...
           Information.subject_info.drug, Information.subject_info.condition];
    fprintf(fileID,'%s \t %s \t %s \t %s \n', sum);
    fclose(fileID); 
    
%     % save CSER
%     T = table(Information.CSER', 'VariableNames', ["CSER Var(1)"]);
%     writetable(T, strcat(path,'CSER.csv'));


end