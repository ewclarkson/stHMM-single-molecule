%% Generate a figure of a "raw data" trajectory

% Assume that D1 > D2

% model parameters
k12 = 8; % association rate % 4, 40, 400
k21 = 8; % dissociation rate
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

% f = figure('Position',[500 200 600 400]);
plot(dataX,dataY, 'LineWidth', 1, 'Marker', 'o','MarkerFaceColor',"#0072BD")
yticks([])
xticks([])
% xlabel('x-coordinate (\mum)')
% ylabel('y-coordinate (\mum)')
% grid on

disp(['nu = ', num2str(max(k12*tau,k21*tau))])

%% colour-code displacements based on w-value

% hold on
% newLen = 30; % to shorten trajectory
newLen = nSteps;
dataX = dataX(1:newLen);
dataY = dataY(1:newLen);

% f = figure('Position',[500 200 600 400]);
% plot(dataX,dataY,'Marker','o','MarkerSize',4,'MarkerFaceColor',"#0072BD")

% cmap = colormap('parula');
cmap = [1 0 0; 0 0 1];
% colormap(cmap)
cimap = [0;1];
cmap2 = interp1(cimap,cmap,linspace(0,1,1E3));
colormap(cmap2)
% wmap = linspace(0, 1, length(cmap));
for i = 1:newLen-1   
    icol = [wMat(i) 0 1-wMat(i)]; % NOTE: update this to cmap2
    % plot(dataX(i:i+1), dataY(i:i+1), 'color', icol,...
    %     'LineWidth', 1, 'Marker', 'o', 'MarkerFaceColor', icol);
    plot(dataX(i:i+1), dataY(i:i+1), 'color', icol,...
        'LineWidth', 2,'Marker', 'o','MarkerFaceColor',icol);
    hold on
end
% xlabel('x-coordinate (\mum)')
% ylabel('y-coordinate (\mum)')
yticks([])
xticks([])
% xlabel(append('\nu = ', num2str(max(k12*tau,k21*tau))))
% grid on
% cbar = colorbar;
% cbar.Label.String = 'w_n';
% c.Label.FontSize = 12;
% cbar.Location = 'east'; % colorbar inside of plot


%% make histograms

% model parameters
k12 = 400; % association rate % 4, 40, 400
k21 = 400; % dissociation rate
D1 = 10; % diffusion constant of state 1
D2 = 1; % diffusion constant of state 2

% data parameters
nSteps = 1E4; % trajectory length
nTraj = 1; % number of trajectories

% experimental parameters
tau = 0.005;
nSub = 100;

[dataX, dataY, stateMat] = util.generatetrajectory_plot(nSteps, tau, nSub, D1, D2, k12, k21);
mask1 = stateMat==1;
wMat = sum(mask1,2);
wMat = wMat/nSub;

% wVec = rand(1,1E4);
wVec = wMat';
wVec = sort(wVec); % used to find w_n values for plotting
[nCounts, edges, bin] = histcounts(wVec); 

nBins = length(edges)-1;
binWidth = edges(2)-edges(1);
% binCentra = edges(1:end-1)+binWidth;
% binCentra = edges(1:end-1)+binWidth/2;
binCentra = linspace(0,1,nBins); % redo?

mybar = bar(binCentra,nCounts,1,'EdgeColor','none');
% xticks(round(binCentra,2))
xticks([0,0.25,0.5,0.75,1])
xlabel('w_n (hidden)')
ylabel('count')

ind = find(diff(bin)~=0); % find "transition points" of bins
binInd = [1 ind]; % double check correctness of indices

mybar.FaceColor = 'flat';
for idx = 1:nBins
    mybar.CData(idx,:) = [1/nBins*idx 0 1-1/nBins*idx]; % correct?
    % mybar.CData(idx,:) = [wVec(binInd(idx)) 0 1-wVec(binInd(idx))]; % better?
end


%%

cmap = [0 0 1; 1 0 0];
cimap = [0;1];
cmap2 = interp1(cimap,cmap,linspace(0,1,1E3));

f1 = openfig('figures/fig1_underlying_nu=0.02.fig');
ax1 = findall(f1,'type','axes');

f2 = openfig('figures/fig1_histogram_nu=0.02.fig');
ax2 = findall(f2,'type','axes');

f3 = openfig('figures/fig1_underlying_nu=0.2.fig');
ax3 = findall(f3,'type','axes');

f4 = openfig('figures/fig1_histogram_nu=0.2.fig');
ax4 = findall(f4,'type','axes');

f5 = openfig('figures/fig1_underlying_nu=2.fig');
ax5 = findall(f5,'type','axes');

f6 = openfig('figures/fig1_histogram_nu=2.fig');
ax6 = findall(f6,'type','axes');

f = figure('Position',[500 100 600 450]);
tl = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

% plot 1
ax1c = copyobj(ax1,tl);
ax1c(1).Layout.Tile = 1;
ax1c.Layout.TileSpan = [1 1];
text(0.01,0.95,{"\bf{a}"},'FontSize',12,'units','normalized');

% plot 2
ax2c = copyobj(ax2,tl);
ax2c(1).Layout.Tile = 2;
ax2c.Layout.TileSpan = [1 1];
text(0.01,0.92,ax2,{"\bf{b}"},'Parent',ax2c,'FontSize',12,'units','normalized');

% plot 3
ax3c = copyobj(ax3,tl);
ax3c(1).Layout.Tile = 3;
ax3c.Layout.TileSpan = [1 1];
text(0.01,0.95,ax3,{"\bf{c}"},'Parent',ax3c,'FontSize',12,'units','normalized');

% plot 4
ax4c = copyobj(ax4,tl);
ax4c(1).Layout.Tile = 4;
ax4c.Layout.TileSpan = [1 1];
text(0.01,0.92,ax4,{"\bf{d}"},'Parent',ax4c,'FontSize',12,'units','normalized');

% plot 5
ax5c = copyobj(ax5,tl);
ax5c(1).Layout.Tile = 5;
ax5c.Layout.TileSpan = [1 1];
text(0.01,0.95,ax5,{"\bf{e}"},'Parent',ax5c,'FontSize',12,'units','normalized');

% plot 6
ax6c = copyobj(ax6,tl);
ax6c(1).Layout.Tile = 6;
ax6c.Layout.TileSpan = [1 1];
text(0.01,0.92,ax6,{"\bf{f}"},'Parent',ax6c,'FontSize',12,'units','normalized');

% % colourbar
% colormap(cmap2)
% cbar2 = colorbar;
% cbar2.Label.String = 'w_n (unobserved)';
% % cbar2.Label.FontSize = 12;
% % cbar2.Location = 'east'; % colorbar inside of plot
% cbar2.Location = 'layout'; % colorbar in its own tile
% cbar2.Layout.Tile = 17;
% cbar2.Layout.TileSpan = [3 1];
