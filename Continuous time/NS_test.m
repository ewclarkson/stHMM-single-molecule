%% test nested sampling for non-uniform priors

mu1 = 7; % mean of state 1
mu2 = 3; % mean of state 2
sigma = 1; % std of both states
k12 = 100; 
k21 = k12;

% nested sampling parameters
nLive = 100; % number of live points
StopRatio = 1E-4; % stop criterion for evidence

% data parameters
nSteps = 300; % trajectory length
nTraj = 10; % number of trajectories

% prior
% priorPars = [4,9; 1,5; 0.001,0.5; 0.001,0.5]; % uniform
% priorPars = {'uniform',4,9; 'uniform',1,5; 'uniform',0.001,0.5; 'uniform',0.001,0.5};
% priorPars = {'lognormal',log(mu1),log(1.5); 'uniform',1,5; 'uniform',0.001,0.5; 'uniform',0.001,0.5}; 
priorPars = {'lognormal',log(mu1),log(1.5); 'lognormal',log(mu2),log(1.5); 'lognormal',log(k12),log(5.5); 'lognormal',log(k21),log(5.5)}; 

% experimental parameters
tau = 5; % sampling time
exPars = {'tau', 5; 'sigma', sigma}; % defines experimental parameters that enter models

% generate data
data = cell(1,nTraj);
for i = 1:nTraj
    data{i} = util.generatesignal(nSteps, tau, 100, 100, mu1, mu2, sigma, sigma, k12, k21);
end

[finalSeq_DTHMM, thetaMLE, logZ_DTHMM] = util.nestedsampling(nLive, StopRatio, priorPars, exPars, data, @util.logl_DTHMM_signal);
[logZ_error_DTHMM, thetaBayes_DTHMM, errorBayes_DTHMM] = util.bayesianestimate(finalSeq_DTHMM,logZ_DTHMM,nLive,0);
thetaBayes_DTHMM

% [finalSeq_CTHMM, thetaMLE, logZ_CTHMM] = util.nestedsampling(nLive, StopRatio, priorPars, exPars, data, @util.logl_CTHMM_signal);
% [logZ_error_CTHMM, thetaBayes_CTHMM, errorBayes_CTHMM] = util.bayesianestimate(finalSeq_CTHMM,logZ_CTHMM,nLive,0);
% thetaBayes_CTHMM


%% plot priors

X = 0.05:0.01:18.95;
mu = 0.5;
sigma = 1;
Y = lognpdf(X,mu,sigma);

plot(log(X),Y)

%% mu1 prior

X = 0.05:0.01:18.95;
mu = log(mu1);
sigma = log(1.5);
Y = lognpdf(X,mu,sigma);

plot(X,Y)

%% mu2 prior

X = 0.05:0.01:18.95;
mu = log(mu2);
sigma = log(1.5);
Y = lognpdf(X,mu,sigma);

plot(X,Y)

%% k12,k21 priors

% X = 1E-6:0.001:1; % units of millisecond
% mu = log(0.1);
X = 1:0.1:600; % units of second
mu = log(200);
sigma = log(3.5);
Y = lognpdf(X,mu,sigma);

plot(X,Y)





