%% Generate a figure of a "raw data" trajectory

% Assume that D1 > D2

% model parameters
k12 = 400; % association rate % 4, 40, 400
k21 = 400; % dissociation rate
D1 = 10; % diffusion constant of state 1
D2 = 1; % diffusion constant of state 2

% data parameters
nSteps = 3; % trajectory length

% experimental parameters
tau = 0.005;
Rmb = 1/6;
sigmaE = 0.010;
exPars = {'tau', tau; 'Rmb', Rmb; 'sigmaE', sigmaE}; % experimental parameters that enter models
nSub = 1E2;

% generate trajectories
% [dataX, dataY] = util.generatetrajectory_pos(nSteps, tau, 100, D1, D2, k12, k21);
[subX, subY, X, Y, stateMat] = ...
    util.generatetrajectory_simplot(nSteps, tau, nSub, D1, D2, k12, k21,sigmaE);

mask1 = stateMat==1;
wMat = mask1'; % to be able to use linear indexing below

cmap = [0 0.4 1; 1 0.1 0];
% colormap(cmap)
cimap = [0;1];
nColors = 1E3;
cmap2 = interp1(cimap,cmap,linspace(0,1,nColors));

for i = 1:length(subX)-1   
    % icol = [wMat(i) 0 1-wMat(i)];
    icol = cmap2(round(1+wMat(i)*(nColors-1)),:);
    plot(subX(i:i+1),subY(i:i+1),'Color',icol,'LineWidth',1.2)
    hold on
end
plot(X,Y,'LineWidth',2,'Color',[0.5 0.5 0.5],...
    'Marker', '+','MarkerSize',12,'MarkerEdgeColor',[0,0,0])

grid on
% xlabel(append('\nu = ', num2str(max(k12*tau,k21*tau))))

% yticks([])
% xticks([])
% xlabel('x-coordinate (\mum)')
% ylabel('y-coordinate (\mum)')
% grid on

disp(['nu = ', num2str(max(k12*tau,k21*tau))])



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
    % mybar.CData(idx,:) = [1/nBins*idx 0 1-1/nBins*idx];
    mybar.CData(idx,:) = [1/(nBins-1)*(idx-1) 0.4-0.3/(nBins-1)*(idx-1) 1-1/(nBins-1)*(idx-1)]; 
    % check correction with cmap above
    % mybar.CData(idx,:) = [wVec(binInd(idx)) 0 1-wVec(binInd(idx))]; % better?
end


%%

f1 = openfig('figures/fig1_sub_nu=0.02_new.fig');
xlabel('\nu = 0.02')
ax1 = findall(f1,'type','axes');

f2 = openfig('figures/fig1_histogram_nu=0.02_new.fig');
ax2 = findall(f2,'type','axes');

f3 = openfig('figures/fig1_sub_nu=0.2_new.fig');
xlabel('\nu = 0.2')
ax3 = findall(f3,'type','axes');

f4 = openfig('figures/fig1_histogram_nu=0.2_new.fig');
ax4 = findall(f4,'type','axes');

f5 = openfig('figures/fig1_sub_nu=2_new.fig');
xlabel('\nu = 2')
ax5 = findall(f5,'type','axes');

f6 = openfig('figures/fig1_histogram_nu=2_new.fig');
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


