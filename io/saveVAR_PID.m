%% saveVAR_PID
% 
% Saves PID atoms and related information from a VAR PID procedure to files.
% 
% Description:
%   This function saves the PID atoms calculated from a VAR PID procedure to CSV files. 
%   It organizes the files into directories based on redundancy functions, subject information, and model parameters. 
%   The function also saves some parameters and subject information to separate text files for reference.
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
%                     - runs      : Number of runs in the VAR PID procedure.
%                     - red_fun   : Array of redundancy functions used (e.g., "MMI", "DEP", "CCS", or "all").
%                 * Red, Syn, UnX, UnY - Cell arrays containing the computed PID atoms for each redundancy function:
%                     - Red  : Redundancy atoms.
%                     - Syn  : Synergy atoms.
%                     - UnX  : Unique Information X atoms.
%                     - UnY  : Unique Information Y atoms.
%   out         - Output directory path where files will be saved.
% 
% Outputs:
%   The function does not return any outputs. It saves the following files to the specified output directory:
%   * data.csv    - CSV files containing the PID atoms (Redundancy, Synergy, Unique Information X, and Unique Information Y) for each redundancy function.
%   * summary.txt - Text file summarizing the analysis parameters (channels, epochs, and number of runs).
%   * ID.txt      - Text file containing the subject's ID, name, drug, and condition information.
% 
% Notes:
%   - The function generates the output path based on the subject's drug condition, name, number of channels, and state.
%   - If the specified directory does not exist, it is created.
%   - For each redundancy function, a separate folder is created, and the corresponding PID atoms are saved in CSV format.
%   - If "all" is specified in `red_fun`, the function creates folders and saves data for all predefined redundancy functions.
%   - The summary and ID information are saved in separate text files for documentation purposes.
% 
% Alberto Liardi, 2024


function [] = saveVAR_PID(Information, out)

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
    
    PIDs = zeros(3,4,length(Information.PIDs));
    for i = 1:length(Information.PIDs), PIDs(:,:,i) = Information.PIDs{i}; end

    for r = 1:length(Information.parameters.red_fun)
        if Information.parameters.red_fun(r) ~= "all"
            final_path = path+Information.parameters.red_fun(r)+'/';
            if not(isfolder(final_path)), mkdir(final_path); end
            
            if Information.parameters.red_fun(r)=="MMI", i=1;
            elseif Information.parameters.red_fun(r)=="DEP", i=2;
            elseif Information.parameters.red_fun(r)=="CCS", i=3;
            end

            T = array2table(squeeze(PIDs(i,:,:))', ...
                  'VariableNames', ["Unique_information_X", ...
                  "Unique_information_Y", "Redundancy", "Synergy"] );
    
            writetable(T, final_path+'data.csv'); % ,'Delimiter','\t'   
        else
            red_funs = ["MMI","DEP","CCS"];
            for j = 1:3
                final_path = path+red_funs(j)+'/';
                if not(isfolder(final_path)), mkdir(final_path); end

                T = array2table(squeeze(PIDs(j,:,:))', ...
                  'VariableNames', ["Unique_information_X", ...
                  "Unique_information_Y", "Redundancy", "Synergy"] );
    
                writetable(T, final_path+'data.csv'); % ,'Delimiter','\t'
            end
            break;
        end
       
    end
    
    % save some more parameters
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