%% test analysing data from the HMJP article

% nested sampling parameters
nLive = 100; % number of live points
StopRatio = 1E-4; % stop criterion for evidence
% priorLimits = [mu1_min,mu1_max;mu2_min,mu2_max;k12_min,k12_max;k21_min,k21_max];
priorPars = {'lognormal',log(7),log(1.5); 'lognormal',log(1),log(1.5); 'lognormal',log(0.2),log(3.5); 'lognormal',log(0.2),log(3.5)}; 

% experimental parameters
% tau = 0.1; % sampling time
sigma = 1;
exPars = {'tau', 0.1; 'sigma', sigma}; % defines experimental parameters that enter models

% load()
data = {observation.obs}; % first load the data from the HMJP code folder

%%
tic
% run nested sampling
[finalSeq, thetaMLE, logZ] = util.nestedsampling(nLive, StopRatio, priorPars, exPars, data, @util.logl_CTHMM_signal);
seqLen = length(finalSeq);

% extract samples
dataEx = zeros(seqLen,4); % all coordinate values
w = zeros(1,seqLen); % corresponding weights
logLVec = zeros(seqLen,1); % likelihoods

for i = 1:seqLen
    
    dataEx(i,:) = finalSeq(i).pos;
    w(i) = finalSeq(i).postWt;
    logLVec(i) = exp(finalSeq(i).logL);
end
data1 = dataEx(:,1); % 1st coordinate
data2 = dataEx(:,2); % 2nd coordinate
data3 = dataEx(:,3); % 3rd coordinate
data4 = dataEx(:,4); % 4th coordinate

% Compute Bayesian estimates and errors
logZ_error = sqrt(w*(logLVec-logZ)/nLive);
thetaBayes = [w*data1,w*data2,w*data3,w*data4]; % posterior means
M2Bayes = [w*(data1.^2),w*(data2.^2),w*(data3.^2),w*(data4.^2)]; % 2nd moment
stdTheta = sqrt(M2Bayes-thetaBayes.^2); % standard deviation of estimates

% thetaMLE % MLE estimate
% thetaBayes % Bayesian estimate
% logZ_error
disp(['estimates with CTHMM: ', num2str(thetaBayes)])
time_stHMM=toc