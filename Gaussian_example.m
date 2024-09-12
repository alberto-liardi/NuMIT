%% GAUSSIAN EXMAPLE USAGE
% Calculate NuMIT-normalised atoms for a Gaussian system with 'red_fun' definition

% choose PID definition
red_fun = "CCS";

if red_fun=="MMI", PID_Gaussian = @PID_MMI_Gaussian;
elseif red_fun=="DEP", PID_Gaussian = @PID_DEP_Gaussian;
elseif red_fun=="CCS", PID_Gaussian = @PID_CCS_Gaussian;
end

% consider a "real" Gaussian system T=A*S+epsilon and its PID atoms
% here we take a as an example a predominantly redundant one with S=2, T=1
A = [0.45, 0.45];
eps = 0.00001;
Sigma_s = [1, 1-eps; 1-eps, 1];
Sigma_eps = 1;

% calculate PID
Sigma_t = A*Sigma_s*A'+Sigma_eps;
Sigma_full = [Sigma_s, (A*Sigma_s)'; A*Sigma_s, Sigma_t];
MI = entropy(Sigma_t)+entropy(Sigma_s)-entropy(Sigma_full);

[UnX, UnY, Red, Syn] = PID_Gaussian(Sigma_full,2,1);

% set the appropriate model 
model = struct('name',"Gauss",'n',100,'S',2,'T',1,'red_fun',red_fun);
% normalise these atoms using NuMIT
[qUnX, qUnY, qRed, qSyn] = NuMIT_PID(UnX, UnY, Red, Syn,MI,model);

% print the results
fprintf("The quantiles obtained are:\nUnique X = %f\nUnique Y = %f," + ...
        "\nRedundancy = %f,\nSynergy = %f\n\n", qUnX, qUnY, qRed, qSyn);



function H = entropy(Sigma)
    arg = det(2*pi*exp(1)*Sigma);    
    H = 0.5*log(arg);
end