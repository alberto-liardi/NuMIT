%% VAR-PID EXAMPLE USAGE 
% Example code of VAR-PID procedure for multivariate time series data:
% fit a VAR model and calculate PID atoms

rng("default"); 

% output directory
out = "Results_PID/";

% data information (subject number/name, number of channels, conditions)
subject = "A";
nb_channels = 4; 
conditions = ["Test A", "Test B"];

for cond = 1:length(conditions)
    subject_info.name = subject;
    subject_info.condition = conditions(cond);
    % sampling frequency of the data, needed for downsampling and lowpassing
    subject_info.fs = 600;

    fprintf("doing subject %s condition %s and channels %d \n", ...
            subject_info.name,subject_info.condition,nb_channels);

    % set the parameters and prepare the data
    parameters = struct('channels', nb_channels, 'runs', 100, 'epochs', 30, ...
                        'mmorder', 1, 'red_fun', ["MMI"]);
    
    % load the data. 
    % here we use synthetic time series generated from a random VAR model 
    % NB: data must be a cell array in which each entry is an epoch of
    % dimensions [time x channels]
    [data, check] = generate_ts(nb_channels);
    while check==1, [data, check] = generate_ts(nb_channels); end

    % run the analysis
    Information = VAR_analysis(data, subject_info, parameters, ["PID"]);

    % saving the results in 'out' folder
    saveVAR_PID(Information, out);
end


%% NuMIT NORMALISATION

% Load the Raw atoms calculated above and normalise them via NuMIT
red_fun = "MMI";
for cond = 1:length(conditions)
    fprintf("Doing subject %s condition %s and %d channels \n", subject, conditions(cond), nb_channels);
    outpath = out+"/"+subject+"_c"+nb_channels+"/"+conditions(cond)+"/"+red_fun;
    
    % load the raw atoms
    [UnX, UnY, Red, Syn, MI] = LoadPID(outpath+"/data.csv");
    
    % set the model and run the null procedure
    model = struct('name',"VAR",'S',nb_channels,'n',100,'p',ones(length(MI)),'red_fun',red_fun);
    [qUnX, qUnY, qRed, qSyn] = NuMIT_PID(UnX, UnY, Red, Syn, MI, model);
    
    % save the quantiles (normalised atoms)
    SaveQuantPID(qUnX, qUnY, qRed, qSyn, MI, outpath, 1);
end


%% PLOTTING RESULTS

%%% Plot results for the Raw atoms %%%
red_fun = "MMI";
for cond=1:2
    % load the raw atoms
    outpath = out+"/"+subject+"_c"+nb_channels+"/"+conditions(cond)+"/"+red_fun;
    [UnX, UnY, Red, Syn, MI] = LoadPID(outpath+"/data.csv");
    pid_atoms(:,:,cond) = [UnX; UnY; Red; Syn];
end
atoms = pid_atoms(:,:,1)-pid_atoms(:,:,2);
% set the labels
labels = ["Unique X", "Unique Y", "Redundancy", "Synergy"];
% plot the differences in a swarmchart plot
ViolinPlot(atoms,labels,conditions(1)+" - "+conditions(2),"Raw atoms");


%%% Plot results for the NuMIT atoms %%%
for cond=1:2
    % load the NuMIT-normalised atoms
    outpath = out+"/"+subject+"_c"+nb_channels+"/"+conditions(cond)+"/"+red_fun;
    [UnX, UnY, Red, Syn, MI] = LoadPID(outpath+"/Quantiles_PID_1.csv");
    % transform cell of struct into a 2D array
    pid_atoms(:,:,cond) = [UnX; UnY; Red; Syn];
end
atoms = pid_atoms(:,:,1)-pid_atoms(:,:,2);
% set the labels
labels = ["Unique X", "Unique Y", "Redundancy", "Synergy"];
% plot the differences in a swarmchart plot
ViolinPlot(atoms,labels,conditions(1)+" - "+conditions(2),"NuMIT atoms");


function [data, check] = generate_ts(nvars)

    check = 0;
    % Generate random VAR coefficients
    A = var_rand(nvars,6,0.95);
    % Generate random residuals covariance matrix.
    V = corr_rand(nvars,0.5);
    % Report information on the generated VAR model and check for errors.
    if var_info(A,V,0).error
        disp('VAR error(s) found - bailing out');
        check = 1; 
        data = nan;
        return;
    end
    
    % Generate multi-trial VAR time for generated VAR coefficients 
    X = var_to_tsdata(A,V,1200,100); % time length 1200 and 100 trials
    data = cell(1,100);
    for d=1:100, data{d} = X(:,:,1)'; end

end