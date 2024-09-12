%% SAMPLE USAGE OF THE DISCRETE NULL MODEL

rng("default");

%%%%%%%%%%%%% FIRST APPROACH %%%%%%%%%%%%% 

% consider a "real" system - e.g. a XOR with equal probabilities - and
% calculates its PID atoms
ps = [0.25, 0.25, 0.25, 0.25];
pe = 0.49;
[UnX, UnY, Red, Syn, MI] = PID_MMI_Discrete(ps,pe,"XOR");

% now construct the null models for the given value of MI
null_pids = Null_model_Discrete(MI);

% NB: MMI assigns zero to either of the unique information, hence the null
% distribution is the combination of the two (to avoid sequences of zeros)
null_pids = [sum(null_pids(1:2,:),1); sum(null_pids(1:2,:),1); null_pids(3:4,:)];

% now take the quantiles of the original PID atoms w.r.t. the null distribution
quant_pids = CompQuantile(null_pids, [UnX, UnY, Red, Syn]');

% print the results
fprintf("The quantiles obtained are:\nUnique X = %f\nUnique Y = %f," + ...
        "\nRedundancy = %f,\nSynergy = %f\n\n", quant_pids);


%%%%%%%%%%%%% SECOND APPROACH %%%%%%%%%%%%% 
% In alternative, call the function NuMIT_PID specifying the model for
% which the quantiles need to be computed

% start from the system above
ps = [0.25, 0.25, 0.25, 0.25];
pe = 0.49;
[UnX, UnY, Red, Syn, MI] = PID_MMI_Discrete(ps,pe,"XOR");

% create 'model' struct with necessary information
model.name = "Discrete";
model.n = []; % use default number for the size of the null distribution

% call NuMIT_PID and obtain the quantiles of the atoms
[qUnX, qUnY, qRed, qSyn] = NuMIT_PID(UnX, UnY, Red, Syn, MI,model);
quant_pids = [qUnX, qUnY, qRed, qSyn];

% print the results
fprintf("The quantiles obtained are:\nUnique X = %f\nUnique Y = %f," + ...
        "\nRedundancy = %f,\nSynergy = %f\n\n", quant_pids);
