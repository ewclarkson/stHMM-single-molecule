%% Generate a figure of a "raw data" signal

% Assume that D1 > D2

% model parameters
k12 = 50; % association rate
k21 = 50; % dissociation rate
mu1 = 7;
mu2 = 3;
sigma = 1;

% data parameters
nSteps = 30; % trajectory length
nTraj = 1; % number of trajectories

% experimental parameters
tau = 0.005;
nSub = 50; % number of sub-time intervals per time interval

% generate trajectories
[dataX,stateMat] = util.generatesignal_plot(nSteps, tau, nSub, nSub, mu1, mu2, sigma, sigma, k12, k21);

f = figure('Position',[500 200 600 200]);
% plot(tau*[1:1:nSteps],dataX,'LineWidth',1.0,'Marker','o','MarkerSize',4,'MarkerFaceColor',"#0072BD")
plot(tau*[1:1:nSteps],dataX,'LineWidth',1.5)
xticks(tau*[1:1:nSteps])
xticklabels({});
% xlabel('time (s)')
ylabel('Signal (a.u.)')
set(gca,'XAxisLocation','top');
xlabel('Time')

%% ground truth state sequence

M2 = stateMat';
M2s = M2(:);
X = [1:length(M2s)];
plot(X,M2s,'Color','red','LineWidth',1.5)
ylabel('State')
xlabel('Time')
yticks([1,2])
xticks([1:nSub:nSteps*nSub]);
xticklabels({});
% set(gca,'XAxisLocation','top');

%% state 1 weights
figure;
mask1 = stateMat==1;
wMat = sum(mask1,2);
wMat = wMat/nSub;

% plot(wMat,'LineWidth',1.5,'Marker','o','MarkerSize',4,'MarkerFaceColor',"#0072BD")
plot(wMat,'LineWidth',1.5)
ylabel('w_n')
xlabel('Time')
xticks(1:length(wMat))
xticklabels({})

