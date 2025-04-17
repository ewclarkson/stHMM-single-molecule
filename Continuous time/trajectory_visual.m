%% Generate a figure of a "raw data" trajectory

% Assume that D1 > D2

% model parameters
k12 = 200; % association rate
k21 = 200; % dissociation rate
D1 = 10; % diffusion constant of state 1
D2 = 1; % diffusion constant of state 2

% data parameters
nSteps = 300; % trajectory length
nTraj = 1; % number of trajectories

% experimental parameters
tau = 0.005;
exPars = {'tau', tau; 'Rmb', 0; 'sigmaE', 0}; % defines experimental parameters that enter models
nSub = 30;

% generate trajectories
% [dataX, dataY] = util.generatetrajectory_pos(nSteps, tau, 100, D1, D2, k12, k21);
[dataX, dataY, stateMat] = util.generatetrajectory_plot(nSteps, tau, nSub, D1, D2, k12, k21);

mask1 = stateMat==1;
wMat = sum(mask1,2);
wMat = wMat/nSub;

f = figure('Position',[500 200 600 200]);
plot(dataX,dataY)
xlabel('x-coordinate (\mum)')
ylabel('y-coordinate (\mum)')
grid on

%% colour-code displacements based on w-value

close all
% shorten the trajectory
newLen = 50;
dataX = dataX(1:newLen);
dataY = dataY(1:newLen);

f = figure('Position',[500 200 600 400]);
% plot(dataX,dataY,'Marker','o','MarkerSize',4,'MarkerFaceColor',"#0072BD")

% cmap = colormap('parula');
cmap = [1 0 0; 0 0 1];
% colormap(cmap)
cimap = [0;1];
cmap2 = interp1(cimap,cmap,linspace(0,1,1E3));
colormap(cmap2)
% wmap = linspace(0, 1, length(cmap));
for i = 1:newLen-1   
    % icol = interp1(wmap, cmap, wMat(i));
    % icol = cmap2(,:);
    icol = [1-wMat(i) 0 wMat(i)]; % NOTE: update this to cmap2
    plot(dataX(i:i+1), dataY(i:i+1), 'color', icol,...
        'LineWidth', 1, 'Marker', 'o', 'MarkerFaceColor', icol);
    % plot(dataX(i:i+1), dataY(i:i+1));
    hold on
end
xlabel('x-coordinate (\mum)')
ylabel('y-coordinate (\mum)')
grid on
% clim([0, 1])
cbar = colorbar;
cbar.Label.String = 'w_n';
c.Label.FontSize = 12;
cbar.Location = 'east'; % colorbar inside of plot
% cbar.Location = 'manual';
% cbar.Position = [0.79 0.1100 0.0356 0.8160];


%% Load and combine multiple figures

f1 = openfig('figures/trajectory_visual_2025-04-15.fig');
% f1.Position = [500 200 600 500];
ax1 = gca;

f2 = openfig('figures/trajectory_main_2025-04-15.fig');
f2.Position = [500 200 600 500];
ax2_all = findall(f2,'type','axes');
lgd = ax2_all(4).Legend;

f3 = figure('Position',[500 100 600 650]);
tl = tiledlayout(8,2,'TileSpacing','compact','Padding','compact');
% tl = tiledlayout(8,6,'TileSpacing','compact','Padding','compact');

ax1c = copyobj(ax1,tl);
ax1c(1).Layout.Tile = 1;
ax1c.Layout.TileSpan = [2 2];
% ax1c.Layout.TileSpan = [2 5];
% nexttile
colormap(cmap2)
cbar2 = colorbar;
cbar2.Label.String = 'w_n';
cbar2.Label.FontSize = 12;
cbar2.Location = 'east'; % colorbar inside of plot
% cbar2.Location = 'layout'; % colorbar in its own tile
% cbar2.Layout.Tile = 6;
% cbar2.Layout.TileSpan = [2 1];

% plot 1
ax2cGroup = copyobj([ax2_all(4),lgd],tl);
ax2c = ax2cGroup(1);
lgdc = ax2cGroup(2);
ax2c(1).Layout.Tile = 4+1;
ax2c.Layout.TileSpan = [3 1];
% ax2c = copyobj(ax2_all(4),tl);
% ax2c(1).Layout.Tile = 4+1;
% ax2c.Layout.TileSpan = [3 1];

% plot 2
ax2c = copyobj(ax2_all(3),tl);
ax2c(1).Layout.Tile = 4+2;
ax2c.Layout.TileSpan = [3 1];

% plot 3
ax2c = copyobj(ax2_all(2),tl);
ax2c(1).Layout.Tile = 11;
ax2c.Layout.TileSpan = [3 1];

% plot 4
ax2c = copyobj(ax2_all(1),tl);
ax2c(1).Layout.Tile = 12;
ax2c.Layout.TileSpan = [3 1];


% cou = 1;
% for idx = 1:4
%     cou = cou+idx;
%     ax2c = copyobj(ax2_all(idx),tl);
%     ax2c(1).Layout.Tile = 4+idx*3;
%     ax2c.Layout.TileSpan = [3 1];
% end



%% test

f1 = openfig('figures/trajectory_visual_2025-02-16.fig');
% f2.Position = [500 200 600 500];
% ax1 = axes(f1);
ax1 = gca;
% lgd = legend(ax2);

f2 = openfig('figures/trajectory_main_10traj_300steps_100livepoints.fig');
f2.Position = [500 200 600 500];
% ax1 = axes(f1);
ax2_all = findall(f2,'type','axes');

f3 = figure('Position',[500 100 600 650]);
tl = tiledlayout(3,2,'TileSpacing','loose','Padding','loose');
% f3 = figure;
% tl = tiledlayout(f3,3,2);

% nexttile([1 2]);
ax1c = copyobj(ax1,tl);
ax1c(1).Layout.Tile = 1;
ax1c.Layout.TileSpan = [1 2];

for idx = 1:4
    ax2c = copyobj(ax2_all(idx),tl);
    ax2c(1).Layout.Tile = 2+idx;
end
% ax1c = copyobj(ax1,tl);
% ax2c = copyobj(ax2,tl);
% ax1c(1).Layout.Tile = 1;
% ax2c(1).Layout.Tile = 2;


