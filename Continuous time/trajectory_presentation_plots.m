%% Generate a figure of a "raw data" trajectory

% Assume that D1 > D2

% model parameters
k12 = 50; % association rate
k21 = 50; % dissociation rate
D1 = 10; % diffusion constant of state 1
D2 = 1; % diffusion constant of state 2

% data parameters
nSteps = 30; % trajectory length
nTraj = 1; % number of trajectories

% experimental parameters
tau = 0.005;
nSub = 50; % number of sub-time intervals per time interval

% generate trajectories
[dataX, dataY, stateMat] = util.generatetrajectory_plot(nSteps, tau, nSub, D1, D2, k12, k21);

f = figure('Position',[500 200 300 300]);
plot(dataX,dataY,'Marker','o','MarkerSize',4,'MarkerFaceColor',"#0072BD")
xlabel('x-coordinate (\mum)')
ylabel('y-coordinate (\mum)')
grid on

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

mask1 = stateMat==1;
wMat = sum(mask1,2);
wMat = wMat/nSub;

% plot(wMat,'LineWidth',1.5,'Marker','o','MarkerSize',4,'MarkerFaceColor',"#0072BD")
plot(wMat,'LineWidth',1.5)
ylabel('w_n')
xlabel('Time')
xticks(1:length(wMat))
xticklabels({})

