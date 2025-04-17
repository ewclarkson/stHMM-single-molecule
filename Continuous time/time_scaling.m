%% Test time difference between stHMM & HMM

% model parameters
k12 = 0.1; % association rate
k21 = 0.1; % dissociation rate
mu1 = 5;
mu2 = 3;
sigma = 1;

nSteps = 1E5; % data length
nDatasets = 1; % number of trajectories

% experimental parameters
tau = 5; % sampling time
exPars = {'tau', 5; 'sigma', sigma}; % defines experimental parameters that enter models

% generate data
data = cell(1,nDatasets);
for i = 1:nDatasets
    data{i} = util.generatesignal(nSteps, tau, 100, 100, mu1, mu2, sigma, sigma, k12, k21);
end

% analyse data
paramArr = [mu1, mu2, k12, k21];
% enslogL_Kinz = util.logl_Kinz_signal(paramArr, exPars, data);
% enslogL_DTHMM = util.logl_DTHMM_signal(paramArr, exPars, data);
tic
enslogL_HMM = util.logl_DTHMM_signal(paramArr, exPars, data);
time_HMM = toc;
tic
enslogL_CTHMM = util.logl_CTHMM_signal(paramArr, exPars, data);
time_stHMM = toc;

% time_HMM
% time_stHMM

quotient = time_stHMM/time_HMM

%% Test N-scaling with trajectory length

% generate data
nDatasets = 1;
nSteps = 1E4;
data = cell(1,nDatasets);
for i = 1:nDatasets
    data{i} = util.generatesignal(nSteps, tau, 100, 100, mu1, mu2, sigma, sigma, k12, k21);
end

nTests = 9;
timeArr = zeros(1,nTests);
for idx = 1:nTests
    nSteps = idx*1000; % trajectory length
    data_loop = {data{1}(1:nSteps)}; % truncate the trajectory
    % tic
    % enslogL_CTHMM = util.logl_CTHMM_signal(paramArr, exPars, data_loop);
    % time_stHMM = toc;
    fhandle = @() util.logl_CTHMM_signal(paramArr, exPars, data_loop);
    time_stHMM = timeit(fhandle);
    timeArr(idx) = time_stHMM;
end

plot(timeArr)


%% Test stHMM time scaling with trajectory length

% model parameters
k12 = 0.1; % association rate
k21 = 0.1; % dissociation rate
mu1 = 5;
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
nLive = 100; % number of live points
StopRatio = 1E-4; % stop criterion for evidence
% priorLimits = [mu1_min,mu1_max;mu2_min,mu2_max;k12_min,k12_max;k21_min,k21_max];
priorPars = {'lognormal',log(mu1),log(1.5); 'lognormal',log(mu2),log(1.5); 'lognormal',log(k12),log(5.5); 'lognormal',log(k21),log(5.5)}; 

% experimental parameters
tau = 5; % sampling time
exPars = {'tau', 5; 'sigma', sigma}; % defines experimental parameters that enter models
nDatasets = 30; % number of trajectories

nTests = 4;
timeArr = zeros(1,nTests);
data = cell(1,nDatasets);
fullLength = 1000;
for i = 1:nDatasets
    data{i} = util.generatesignal(fullLength, tau, 100, 100, mu1, mu2, sigma, sigma, k12, k21);
end

for idx = 1:nTests % do for trajectory lengths
    %tempArr = zeros(1,3);
    % data parameters
    nSteps = 50 + idx*50; % trajectory length    
    tic
    for idy = 1:nDatasets
        data_loop = {data{idy}(1:nSteps)}; % truncate the trajectory
        % run nested sampling
        [finalSeq, thetaMLE, logZ] = util.nestedsampling(nLive, StopRatio, priorPars, exPars, data_loop, @util.logl_CTHMM_signal);
        seqLen = length(finalSeq);
    end
    runTime = toc;
    %tempArr(idy) = runTime;
    timeArr(idx) = runTime;
    %timeArr(idx) = mean(tempArr);
    idx
end

timeArr        
plot(timeArr);