%% test generating and analysing a single-molecule signal trajectory

% model parameters
k12 = 0.1; % association rate
k21 = 0.1; % dissociation rate
mu1 = 7;
mu2 = 3;
sigma = 1;

% priors parameters
mu1_max = 20; % known upper limit of D1
mu1_min = 2; % known lower limit of D1
mu2_max = 9; % known upper limit of D2
mu2_min = 1; % known lower limit of D2
k12_min = 1E-3;
k12_max = 5*k12;
k21_min = 1E-3;
k21_max = 5*k21;

% nested sampling parameters
nLive = 300; % number of live points
StopRatio = 1E-4; % stop criterion for evidence
% priorLimits = [mu1_min,mu1_max;mu2_min,mu2_max;k12_min,k12_max;k21_min,k21_max];
priorPars = {'lognormal',log(mu1),log(1.5); 'lognormal',log(mu2),log(1.5); 'lognormal',log(0.2),log(3.5); 'lognormal',log(0.2),log(3.5)}; 

% data parameters
nSteps = 3E2; % data length
nDatasets = 10; % number of trajectories

% experimental parameters
tau = 5; % sampling time
exPars = {'tau', 5; 'sigma', sigma}; % defines experimental parameters that enter models

% generate data
data = cell(1,nDatasets);
for i = 1:nDatasets
    data{i} = util.generatesignal(nSteps, tau, 100, 100, mu1, mu2, sigma, sigma, k12, k21);
end

% analyse data
% paramArr = [mu1, mu2, k12, k21];
% enslogL_Kinz = util.logl_Kinz_signal(paramArr, exPars, data);
% enslogL_DTHMM = util.logl_DTHMM_signal(paramArr, exPars, data);

disp(['nu = ', num2str(max(k12*tau,k21*tau))])


%% ------------------ DTHMM likelihood ---------------------
tic
% run nested sampling
[finalSeq, thetaMLE, logZ] = util.nestedsampling(nLive, StopRatio, priorPars, exPars, data, @util.logl_DTHMM_signal);
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
disp(['estimates with DTHMM: ', num2str(thetaBayes)])
time_HMM = toc
%% ------------------ CTHMM likelihood ---------------------
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
%% ------------------ Kinz likelihood ---------------------

% run nested sampling
[finalSeq, thetaMLE, logZ] = util.nestedsampling(nLive, StopRatio, priorPars, exPars, data, @util.logl_Kinz_signal);
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
disp(['estimates with Kinz: ', num2str(thetaBayes)])


%% Plot marginal distributions

nBins = round(sqrt(seqLen)); % number of bins

figure;
[xbar,binCounts] = util.bincounts(nBins, thetaBayes(1), stdTheta(1), data1, w);
bar(xbar,binCounts,1) % weighted histogram
xlabel('\mu_1')
ylabel('P(\mu_1)')

figure;
[xbar,binCounts] = util.bincounts(nBins, thetaBayes(2), stdTheta(2), data2, w);
bar(xbar,binCounts,1) % weighted histogram
xlabel('\mu_2')
ylabel('P(\mu_2)')

figure;
[xbar,binCounts] = util.bincounts(nBins, thetaBayes(3), stdTheta(3), data3, w);
bar(xbar,binCounts,1) % weighted histogram
xlabel('k_{12}')
ylabel('P(k_{12})')

figure;
[xbar,binCounts] = util.bincounts(nBins, thetaBayes(4), stdTheta(4), data4, w);
bar(xbar,binCounts,1) % weighted histogram
xlabel('k_{21}')
ylabel('P(k_{21})')