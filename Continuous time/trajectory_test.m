%% Compare discrete to continous-time HMM for parameter inference

% Assume that D1 > D2

% model parameters
k12 = 40; % association rate
k21 = 40; % dissociation rate
D1 = 10; % diffusion constant of state 1
D2 = 1; % diffusion constant of state 2

% data parameters
nSteps = 300; % trajectory length
nTraj = 10; % number of trajectories

% experimental parameters
sigmaE = 0.010; %0.015;
Rmb = 1/6; %1/6; % =0: no motion blur, =1/6: motion blur with continuous illumination
tau = 0.005;
exPars = {'tau', tau; 'Rmb', Rmb; 'sigmaE', sigmaE}; % defines experimental parameters that enter models

% generate trajectories
data = cell(1,nTraj);
for i = 1:nTraj
    % data{i} = util.generatetrajectory(nSteps, tau, 100, D1, D2, k12, k21);
    data{i} =  util.generatetrajectory_noisy(nSteps, tau, 1000, 1000, D1, D2, k12, k21, sigmaE);
end

% logL_1state = util.logl_1state([0.5*(D1+D2)],exPars,data)
% logL_DTHMM = util.logl_DTHMM([D1,D2,k12,k21],exPars,data)
% logL_CTHMM = util.logl_CTHMM([D1,D2,k12,k21],exPars,data)

% uniform priors parameters
% D1_max = 20; % known upper limit of D1
% D1_min = 2; % known lower limit of D1
% D2_max = 5; % known upper limit of D2
% D2_min = 0.01; % known lower limit of D2
% k12_min = 10;
% k12_max = 3*k12;
% k21_min = 10;
% k21_max = 3*k21;

% algorithm parameters
nLive = 200; % number of live points
StopRatio = 1E-4; % stop criterion for evidence
% priorPars = {'uniform',D1_min,D1_max; 'uniform',D2_min,D2_max; 'uniform',k12_min,k12_max; 'uniform',k21_min,k21_max};
priorPars = {'lognormal',log(D1),log(1.5); 'lognormal',log(D2),log(1.5);...
             'lognormal',log(200),log(3.5); 'lognormal', log(200),log(3.5)}; 

disp(['nu = ', num2str(max(k12*tau,k21*tau))])

%% ------------------ 1-state diffusion ---------------------

priorPars1state = {'lognormal',log(0.5*(D1+D2)),log(1.5)};

% run nested sampling - 1-state
[finalSeq, thetaMLE, logZ] = util.nestedsampling_1state(nLive, StopRatio, priorPars, exPars, data, @util.logl_1state);
seqLen = length(finalSeq);

% extract samples
dataEx = zeros(seqLen,1); % all coordinate values
w = zeros(1,seqLen); % corresponding weights
logLVec = zeros(seqLen,1); % likelihoods

for i = 1:seqLen
    
    dataEx(i,:) = finalSeq(i).pos;
    w(i) = finalSeq(i).postWt;
    logLVec(i) = exp(finalSeq(i).logL);
end
data1 = dataEx(:,1); % 1st coordinate

% Compute Bayesian estimates and errors
logZ_error = sqrt(w*(logLVec-logZ)/nLive);
thetaBayes = [w*data1]; % posterior means
M2Bayes = [w*(data1.^2)]; % 2nd moment
stdTheta = sqrt(M2Bayes-thetaBayes.^2); % standard deviation of estimates

% thetaMLE % MLE estimate
% thetaBayes % Bayesian estimate
% logZ_error
disp(['estimate with 1-state model is: ', num2str(thetaBayes)])

%% ------------------ DT-HMM ---------------------

% run nested sampling - Das
[finalSeq, thetaMLE, logZ] = util.nestedsampling(nLive, StopRatio, priorPars, exPars, data, @util.logl_DTHMM);
seqLen = length(finalSeq);

% extract samples
dataEx = zeros(seqLen,4); % all coordinate values
w = zeros(1,seqLen); % corresponding weights
logLVec = zeros(seqLen,1); % likelihoods

for i = 1:seqLen
    
    dataEx(i,:) = finalSeq(i).pos;
    w(i) = finalSeq(i).postWt;
    logLVec(i) = finalSeq(i).logL;
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
disp(['estimates with DT-HMM are: ', num2str(thetaBayes)])
disp(['evidence with DT-HMM is: ', num2str(logZ), ' with estimated std: ', num2str(logZ_error)])

%% ------------- CT-HMM --------------------

tic
% run nested sampling - CT-HMM
[finalSeq,thetaMLE,logZ] = util.nestedsampling(nLive, StopRatio, priorPars, exPars, data, @util.logl_CTHMM);
seqLen = length(finalSeq);

% extract data
dataEx = zeros(seqLen,4); % all coordinate values
w = zeros(1,seqLen); % corresponding weights

for i = 1:seqLen

    dataEx(i,:) = finalSeq(i).pos;
    w(i) = finalSeq(i).postWt;
end
data1 = dataEx(:,1); % 1st coordinate
data2 = dataEx(:,2); % 2nd coordinate
data3 = dataEx(:,3); % 3rd coordinate
data4 = dataEx(:,4); % 4th coordinate

% Compute Bayesian estimates and errors
thetaBayes = [w*data1,w*data2,w*data3,w*data4]; % posterior means
M2Bayes = [w*(data1.^2),w*(data2.^2),w*(data3.^2),w*(data4.^2)]; % 2nd moment
stdTheta = sqrt(M2Bayes-thetaBayes.^2); % standard deviation of estimates

% thetaMLE % MLE estimate
disp(['estimates with CT-HMM are: ', num2str(thetaBayes)])
disp(['evidence with CT-HMM is: ', num2str(logZ), ' with estimated std: ', num2str(logZ_error)])
toc


%% Plot marginal distributions

nBins = round(sqrt(seqLen)); % number of bins

figure;
[xbar,binCounts] = util.bincounts(nBins, thetaBayes(1), stdTheta(1), data1, w);
bar(xbar,binCounts,1) % weighted histogram
xlabel('D_1')
ylabel('P(D_1)')

figure;
[xbar,binCounts] = util.bincounts(nBins, thetaBayes(2), stdTheta(2), data2, w);
bar(xbar,binCounts,1) % weighted histogram
xlabel('D_2')
ylabel('P(D_2)')

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

% alpha = 0.025;
% [lb,ub] = util.ci(xbar,binCounts,alpha,alpha)