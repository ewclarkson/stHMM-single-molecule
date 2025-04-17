%% Generate a figure of a "raw data" signal

% Assume that D1 > D2

% model parameters
k12 = 200; % association rate
k21 = 200; % dissociation rate
mu1 = 7;
mu2 = 3;
sigma = 1;

% data parameters
nSteps = 300; % trajectory length
nTraj = 1; % number of trajectories

% experimental parameters
tau = 0.005;
nSub = 30;

% generate trajectories
% dataX = util.generatesignal(nSteps, tau, 100, 100, mu1, mu2, sigma, sigma, k12, k21);
[dataX,stateMat] = util.generatesignal_plot(nSteps, tau, nSub, nSub, mu1, mu2, sigma, sigma, k12, k21);

% calculate weights
mask1 = stateMat==1;
wMat = sum(mask1,2);
wMat = wMat/nSub;

f = figure('Position',[500 200 600 200]);
plot(tau*[0:1:nSteps-1],dataX)
xticks(50*tau*[0:1:nSteps-1])
xlabel('time (s)')
ylabel('signal (a.u.)')
% grid on

%% generate a "zoom-in" and plot time intervals

close all
f = figure('Position',[500 300 600 300]);
% tlZoom = tiledlayout(10,1,'TileSpacing','compact','Padding','compact');
% nexttile([6,1])
% figure

nZoom = 30;
% plot([1:1:nZoom],dataX(1:nZoom),'LineWidth',1.0,...
%     'Marker','o','MarkerSize',3,'MarkerFaceColor',"#0072BD")
% xlabel('time (s)')
% ylabel('signal (a.u.)')
% ------- testing colours ---------
cmap = [0 0 1; 1 0 0];
% colormap(cmap)
cimap = [0;1];
cmap2 = interp1(cimap,cmap,linspace(0,1,1E3));
colormap(cmap2)
for i = 1:nZoom-1
    plot([i,i+1],dataX(i:i+1),'LineWidth',1.0, 'Color', 'black',...
        'Marker','o','MarkerSize',10,'MarkerFaceColor',[wMat(i) 0 1-wMat(i)])
    hold on
end
cbar = colorbar;
cbar.Label.String = 'w_n';
c.Label.FontSize = 12;
xlabel('time (s)')
ylabel('signal (a.u.)')
% --------------------------------

xlVec = [1:1:nZoom];
lblCell = cell(1,nZoom);
wCell = cell(1,nZoom);
for i = 1:nZoom
    lblCell{i} = append(num2str(i),'\tau');
    % wCell{i} = append('w = ',num2str(0.3));
    wCell{i} = append('w = ', num2str(round(wMat(i),1)));
end
xticks(xlVec);
xticklabels(lblCell)
% xline(xlVec,'--', wCell, 'LabelHorizontalAlignment', 'center',...
      % 'LabelVerticalAlignment' ,'top');
% xline(xlVec,'--');
set(gca, 'XGrid','on', 'GridLineStyle','--',...
    'GridLinewidth',1.0, 'GridColor', 'black', 'GridColorMode','manual')

%% plot state sequence of zoom-in
% nexttile([4,1])
figure
% M2 = stateMat(1:nZoom)';
M2 = stateMat';
M2s = M2(:);
M2s = M2s(1:nSub*nZoom);
% plot(tau*[0:1:nZoom-1],M2s,'Color','red','LineWidth',1.5)
plot([1:1:nSub*nZoom],M2s,'Color','red','LineWidth',1.5)
% set(gca, 'YDir','reverse')
ylabel('state')
xlabel('time (s)')
yticks([1,2])
xticks([0:nSub:nZoom*nSub]+1);
% xticklabels({});
% xticks(xlVec);
% xticks(tau*[1:nSub:nZoom*nSub])
xlVec2 = [1:nSub:nZoom*nSub+1];
lblCell2 = cell(1,nZoom+1);
for i = 1:nZoom+1
    lblCell2{i} = append(num2str(i-1),'\tau');
end
xticklabels(lblCell2)
xline(xlVec2,'--')
% set(gca, 'XGrid','on', 'GridLineStyle','--',...
%     'GridLinewidth',1.0, 'GridColor', 'black', 'GridColorMode','manual')

%% Load and combine multiple figures

f1 = openfig('figures/signal_visual_2025-04-11.fig');
% f1.Position = [500 200 600 500];
ax1 = gca;
% lgd = legend(ax2);
fcbar = ax1.Colorbar;

f2 = openfig('figures/signal_main_10traj_200steps_100livepoints.fig');
f2.Position = [500 200 600 500];
ax2_all = findall(f2,'type','axes');
% lgd = legend(ax2_all(4));
lgd = ax2_all(4).Legend;

f3 = figure('Position',[500 100 600 650]);
tl = tiledlayout(8,2,'TileSpacing','compact','Padding','compact');

ax1c = copyobj(ax1,tl);
% ax1c = copyobj([ax1,fcbar],tl);
ax1c(1).Layout.Tile = 1;
ax1c.Layout.TileSpan = [2 2];
colormap(cmap2)
cbar2 = colorbar;
cbar2.Label.String = 'w_n';
c.Label.FontSize = 12;
cbar2.Location = 'east'; % colorbar inside of plot

% plot 1
% ax2c = copyobj(ax2_all(4),tl);
ax2cGroup = copyobj([ax2_all(4),lgd],tl);
ax2c = ax2cGroup(1);
lgdc = ax2cGroup(2);
ax2c(1).Layout.Tile = 4+1;
ax2c.Layout.TileSpan = [3 1];
% ax2c.Legend = lgd;
% legend(lgd)
% legend(ax2c,{'BIASD','HMM','stHMM'},'Location','northwest')

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


