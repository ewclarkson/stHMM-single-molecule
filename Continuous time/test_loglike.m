%% Compare discrete to continous-time HMM for parameter inference

% Assume that D1 > D2

% model parameters
nMul = 100;
k12 = nMul*0.2; % association rate
k21 = nMul*0.2; % dissociation rate
D1 = 10; % diffusion constant of state 1
D2 = 1; % diffusion constant of state 2
paramArr = [D1,D2,k12,k21]; % input into likelihood function

% data parameters
nSteps = 1E2; % trajectory length
nTraj = 10; % number of trajectories

% experimental parameters
exPars = {'tau', 5; 'Rmb', 0; 'sigmaE', 0}; % defines experimental parameters that enter models
tau = 5;

% generate trajectories
data = cell(1,nTraj);
for i = 1:nTraj
    % data{i} = util.generatetrajectory(k12,k21,D1,D2,tau,nSteps);
    data{i} = util.noisyBrownian2state2D(nSteps, tau, 100, 100, D1, D2, k12, k21, 0);
end

% priors parameters
% D1_max = 20; % known upper limit of D1
% D1_min = 2; % known lower limit of D1
% D2_max = 5; % known upper limit of D2
% D2_min = 0.01; % known lower limit of D2
% k12_min = 1E-3;
% k12_max = 3*k12;
% k21_min = 1E-3;
% k21_max = 3*k21;

% algorithm parameters
nLive = 100; % number of live points
StopRatio = 1E-4; % stop criterion for evidence
% priorLimits = [D1_min,D1_max;D2_min,D2_max;k12_min,k12_max;k21_min,k21_max];
priorPars = {'lognormal',log(D1),log(1.5); 'lognormal',log(D2),log(1.5); 'lognormal',log(k12),log(5.5); 'lognormal',log(k21),log(5.5)}; 

disp(['nu = ', num2str(max(k12*tau,k21*tau))])

loglikeCTHMM = util.logl_CTHMM(paramArr, exPars, data);
loglikeKinz = util.logl_Kinz(paramArr, exPars, data);

loglikeCTHMM-loglikeKinz









