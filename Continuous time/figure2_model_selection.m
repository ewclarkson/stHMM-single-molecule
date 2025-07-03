%% Generate a figure of a "raw data" trajectory

% Assume that D1 > D2

% model parameters
k12 = 4; % association rate % 4, 40, 400
k21 = 4; % dissociation rate
D1 = 10; % diffusion constant of state 1
D2 = 1; % diffusion constant of state 2

% data parameters
nSteps = 30; % trajectory length
nTraj = 1; % number of trajectories

% experimental parameters
tau = 0.005;
exPars = {'tau', tau; 'Rmb', 0; 'sigmaE', 0}; % defines experimental parameters that enter models
nSub = 100;

% generate trajectories
% [dataX, dataY] = util.generatetrajectory_pos(nSteps, tau, 100, D1, D2, k12, k21);
[dataX, dataY, stateMat] = util.generatetrajectory_plot(nSteps, tau, nSub, D1, D2, k12, k21);

mask1 = stateMat==1;
wMat = sum(mask1,2);
wMat = wMat/nSub;

f = figure('Position',[500 200 600 400]);
plot(dataX,dataY,'Marker','o','MarkerFaceColor',"#0072BD")
xlabel('x-coordinate (\mum)')
ylabel('y-coordinate (\mum)')
grid on

disp(['nu = ', num2str(max(k12*tau,k21*tau))])

%% colour-code displacements based on w-value

close all
% shorten the trajectory
newLen = 30;
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
    icol = [wMat(i) 0 1-wMat(i)]; % NOTE: update this to cmap2
    plot(dataX(i:i+1), dataY(i:i+1), 'color', icol,...
        'LineWidth', 1, 'Marker', 'o', 'MarkerFaceColor', icol);
    % plot(dataX(i:i+1), dataY(i:i+1));
    hold on
end
% xlabel('x-coordinate (\mum)')
% ylabel('y-coordinate (\mum)')
% yticks([])
% xticks([])
xlabel(append('\nu = ', num2str(max(k12*tau,k21*tau))))
grid on

% cbar = colorbar;
% cbar.Label.String = 'w_n (unobserved)';
% c.Label.FontSize = 12;
% cbar.Location = 'east'; % colorbar inside of plot


%% Load and combine multiple figures

cmap = [0 0 1; 1 0 0];
cimap = [0;1];
cmap2 = interp1(cimap,cmap,linspace(0,1,1E3));

f1 = openfig('figures/visual_trajectory_nu=0.04_2025-05-30.fig');
ax1 = findall(f1,'type','axes');
% ax1 = gca;

f2 = openfig('figures/visual_trajectory_nu=0.5_2025-05-30.fig');
ax2 = findall(f2,'type','axes');

f3 = openfig('figures/visual_trajectory_nu=2.0_2025-05-30.fig');
ax3 = findall(f3,'type','axes');

f4 = openfig('figures/model_probabilities_2025_05_29.fig');
% ax4 = gca;
ax4 = findall(f4,'type','axes');
% f2.Position = [500 200 600 500];

f = figure('Position',[500 100 600 450]);
tl = tiledlayout(10,19,'TileSpacing','compact','Padding','compact');

% plot 1
ax1c = copyobj(ax1,tl);
ax1c(1).Layout.Tile = 1;
ax1c.Layout.TileSpan = [3 5];

% plot 2
ax2c = copyobj(ax2,tl);
ax2c(1).Layout.Tile = 6;
ax2c.Layout.TileSpan = [3 5];

% plot 3
ax3c = copyobj(ax3,tl);
ax3c(1).Layout.Tile = 11;
ax3c.Layout.TileSpan = [3 5];

% colourbar
colormap(cmap2)
cbar2 = colorbar;
cbar2.Label.String = 'w_n (unobserved)';
% cbar2.Label.FontSize = 12;
% cbar2.Location = 'east'; % colorbar inside of plot
cbar2.Location = 'layout'; % colorbar in its own tile
cbar2.Layout.Tile = 17;
cbar2.Layout.TileSpan = [3 1];

% plot 4
ax4c = copyobj(ax4(1),tl);
ax4c(1).Layout.Tile = 58;
ax4c.Layout.TileSpan = [7 19];


%% test

cmap = [0 0 1; 1 0 0];
cimap = [0;1];
cmap2 = interp1(cimap,cmap,linspace(0,1,1E3));

f1 = openfig('figures/visual_trajectory_nu=0.04_2025-05-30.fig');
ax1 = findall(f1,'type','axes');
% ax1 = gca;

f2 = openfig('figures/visual_trajectory_nu=0.5_2025-05-30.fig');
ax2 = findall(f2,'type','axes');

f3 = openfig('figures/visual_trajectory_nu=2.0_2025-05-30.fig');
ax3 = findall(f3,'type','axes');

f4 = openfig('figures/model_probabilities_2025_05_29.fig');
% ax4 = gca;
ax4 = findall(f4,'type','axes');
% f2.Position = [500 200 600 500];

f = figure('Position',[500 100 600 450]);
tl = tiledlayout(10,7,'TileSpacing','compact','Padding','compact');

% plot 1
ax1c = copyobj(ax1,tl);
ax1c(1).Layout.Tile = 1;
ax1c.Layout.TileSpan = [3 2];

% plot 2
ax2c = copyobj(ax2,tl);
ax2c(1).Layout.Tile = 3;
ax2c.Layout.TileSpan = [3 2];

% plot 3
ax3c = copyobj(ax3,tl);
ax3c(1).Layout.Tile = 5;
ax3c.Layout.TileSpan = [3 2];

% colourbar
colormap(cmap2)
cbar2 = colorbar;
cbar2.Label.String = 'w_n (unobserved)';
% cbar2.Label.FontSize = 12;
% cbar2.Location = 'east'; % colorbar inside of plot
cbar2.Location = 'layout'; % colorbar in its own tile
cbar2.Layout.Tile = 7;
cbar2.Layout.TileSpan = [3 1];

% plot 4
ax4c = copyobj(ax4(1),tl);
ax4c(1).Layout.Tile = 22;
ax4c.Layout.TileSpan = [7 7];

%% test

cmap = [0 0 1; 1 0 0];
cimap = [0;1];
cmap2 = interp1(cimap,cmap,linspace(0,1,1E3));

f1 = openfig('figures/visual_trajectory_nu=0.1.fig');
ax1 = gca;

f2 = openfig('figures/visual_trajectory_nu=0.5.fig');
ax2 = gca;

f3 = openfig('figures/visual_trajectory_nu=2.0.fig');
ax3 = gca;

f4 = openfig('figures/model_probabilities.fig');
ax4 = gca;
% f2.Position = [500 200 600 500];

f = figure('Position',[500 100 600 450]);
tl = tiledlayout(10,3,'TileSpacing','compact','Padding','compact');

% colormap(cmap2)
% cbar2 = colorbar;
% cbar2.Label.String = 'w_n (unobserved)';
% cbar2.Location = 'east'; % colorbar inside of plot

% plot 1
ax1c = copyobj(ax1,tl);
ax1c(1).Layout.Tile = 1;
ax1c.Layout.TileSpan = [3 1];

% plot 2
ax2c = copyobj(ax2,tl);
ax2c(1).Layout.Tile = 2;
ax2c.Layout.TileSpan = [3 1];

% plot 3
ax3c = copyobj(ax3,tl);
ax3c(1).Layout.Tile = 3;
ax3c.Layout.TileSpan = [3 1];

% plot 4
ax4c = copyobj(ax4(1),tl);
ax4c(1).Layout.Tile = 10;
ax4c.Layout.TileSpan = [7 3];


