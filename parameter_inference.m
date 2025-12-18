%% Compare discrete to continous-time HMM for parameter inference

% assume D1 > D2
tic
nRep = 7; % try this many values of k12,k21

% model parameters
D1 = 10; % diffusion constant of state 1
D2 = 1; % diffusion constant of state 2
D1vec = linspace(D1,D1,nRep);
D2vec = linspace(D2,D2,nRep);
k12_min_gen = 4;
k12_max_gen = 400;
k21_min_gen = 4;
k21_max_gen = 400;
k12vec = logspace(log10(k12_min_gen),log10(k12_max_gen),nRep);
k21vec = logspace(log10(k21_min_gen),log10(k21_max_gen),nRep);

% set priors
% priorPars = {'lognormal',log(D1),log(1.5); 'lognormal',log(D2),log(1.5); 'lognormal',log(k12vec(i)),log(5.5); 'lognormal',log(k21vec(i)),log(5.5)};
% priorPars = {'uniform',5,15; 'uniform',0.2,2; 'uniform',1,600; 'uniform',1,600};
priorPars = {'lognormal',log(D1),log(1.5); 'lognormal',log(D2),log(1.5); 'lognormal',log(200),log(3.5); 'lognormal',log(200),log(3.5)}; 
% NOTE: k12,k21-priors should be set equal for all data sets

% nested sampling parameters
nLive = 200; % number of live points
StopRatio = 1E-4; % stop criterion for evidence

% data parameters
nSteps = 300; % trajectory length
nTraj = 10; % number of trajectories

% experimental parameters
sigmaE = 0.010; %0.015; % localisation error std
% M = 100; % number of sub-time intervals for trajectory generation
% mC = 1; % motion blur related. 1: no motion blur
Rmb = 1/6; %1/6; % motion blur coefficient. 1/6: continous shutter mode. 0: no motion blur.
tau = 0.005; % sampling time
exPars = {'tau', tau; 'Rmb', Rmb; 'sigmaE', sigmaE}; % defines experimental parameters that enter models
% NOTE: set Rmb to 1/6 if motion blur is included, otherwise to 0!

% nRows = 5; % number of rows, depending on how many values k12-k21 combinations to test

logZcell_DTHMM = cell(1,nRep); % hold estimates of log-evidence
logZcell_Kinz = cell(1,nRep); % hold estimates of log-evidence
logZcell_CTHMM = cell(1,nRep); % hold estimates of log-evidence
theta_DTHMM = cell(1,nRep); % to hold errors of parameter estimates
theta_Kinz = cell(1,nRep); % to hold parameter estimates
theta_CTHMM = cell(1,nRep); % to hold parameter estimates
error_DTHMM = cell(1,nRep); % to hold errors of parameter estimates
error_Kinz = cell(1,nRep); % to hold errors of parameter estimates
error_CTHMM = cell(1,nRep); % to hold errors of parameter estimates

alpha = 0.025; % follows from chosen confidence level
lbounds_DTHMM = zeros(nRep,4); % lower confidence bounds
ubounds_DTHMM = zeros(nRep,4); % upper confidence bounds
lbounds_Kinz = zeros(nRep,4); % lower confidence bounds
ubounds_Kinz = zeros(nRep,4); % upper confidence bounds
lbounds_CTHMM = zeros(nRep,4); % lower confidence bounds
ubounds_CTHMM = zeros(nRep,4); % upper confidence bounds

counter = 0;
for i = 1:nRep
  
    % generate trajectories NOTE: increase the number of sub-time
    % intervals until "convergence"
    data = cell(1,nTraj);
    for k = 1:nTraj
        % data{k} = util.generatetrajectory(nSteps, tau, 100, D1, D2, k12vec(i), k21vec(i));
        data{k} =  util.generatetrajectory_noisy(nSteps, tau, 1000, 1000, D1, D2, k12vec(i), k21vec(i), sigmaE);
    end
    
    % nested sampling with DT-HMM
    [finalSeq_DTHMM, thetaMLE_DTHMM, logZ_DTHMM] = util.nestedsampling(nLive, StopRatio, priorPars, exPars, data, @util.logl_DTHMM);
    
    % run nested sampling with 1-state diffusion
    % [finalSeq_Kinz,thetaMLE_Kinz,logZ_Kinz] = util.nestedsampling(nLive, StopRatio, priorPars, exPars, data, @util.logl_Kinz);

    % run nested sampling with CT-HMM
    [finalSeq_CTHMM,thetaMLE_CTHMM,logZ_CTHMM] = util.nestedsampling(nLive, StopRatio, priorPars, exPars, data, @util.logl_CTHMM);
   
    % Compute Bayesian parameter estimates
    [logZ_error_DTHMM, thetaBayes_DTHMM, errorBayes_DTHMM] = util.bayesianestimate(finalSeq_DTHMM,logZ_DTHMM,nLive,0);
    % [logZ_error_Kinz, thetaBayes_Kinz, errorBayes_Kinz] = util.bayesianestimate(finalSeq_Kinz,logZ_Kinz,nLive,0);
    [logZ_error_CTHMM, thetaBayes_CTHMM, errorBayes_CTHMM] = util.bayesianestimate(finalSeq_CTHMM,logZ_CTHMM,nLive,0);
    
    % store NS output for plotting
    theta_DTHMM{i} = thetaBayes_DTHMM;
    error_DTHMM{i} = errorBayes_DTHMM;
    % theta_Kinz{i} = thetaBayes_Kinz;
    % error_Kinz{i} = errorBayes_Kinz;
    theta_CTHMM{i} = thetaBayes_CTHMM;
    error_CTHMM{i} = errorBayes_CTHMM;

    logZcell_DTHMM{i} = logZ_DTHMM;
    % logZcell_Kinz{i} = logZ_Kinz;
    logZcell_CTHMM{i} = logZ_CTHMM;
    
    [nBins_DTHMM,~,NSpos_DTHMM,w_DTHMM] = util.NSaxuiliaryPlot(finalSeq_DTHMM,logZ_DTHMM,nLive);
    % [nBins_Kinz,~,NSpos_Kinz,w_Kinz] = util.NSaxuiliaryPlot(finalSeq_Kinz,logZ_Kinz,nLive);
    [nBins_CTHMM,~,NSpos_CTHMM,w_CTHMM] = util.NSaxuiliaryPlot(finalSeq_CTHMM,logZ_CTHMM,nLive);
    
    for idx = 1:4 % NOTE: put this loop inside confidenceinterval.m
        [xbar_DTHMM,binCounts_DTHMM] = util.bincounts(nBins_DTHMM, thetaBayes_DTHMM(idx), errorBayes_DTHMM(idx), NSpos_DTHMM(:,idx), w_DTHMM);
        % [xbar_Kinz,binCounts_Kinz] = util.bincounts(nBins_Kinz, thetaBayes_Kinz(idx), errorBayes_Kinz(idx), NSpos_Kinz(:,idx), w_Kinz);
        [xbar_CTHMM,binCounts_CTHMM] = util.bincounts(nBins_CTHMM, thetaBayes_CTHMM(idx), errorBayes_CTHMM(idx), NSpos_CTHMM(:,idx), w_CTHMM);
        [lb_DTHMM,ub_DTHMM] = util.ci(xbar_DTHMM,binCounts_DTHMM,alpha,alpha);
        % [lb_Kinz,ub_Kinz] = util.ci(xbar_Kinz,binCounts_Kinz,alpha,alpha);
        [lb_CTHMM,ub_CTHMM] = util.ci(xbar_CTHMM,binCounts_CTHMM,alpha,alpha);
        lbounds_DTHMM(i,idx) = lb_DTHMM;
        ubounds_DTHMM(i,idx) = ub_DTHMM;
        % lbounds_Kinz(i,idx) = lb_Kinz;
        % ubounds_Kinz(i,idx) = ub_Kinz;
        lbounds_CTHMM(i,idx) = lb_CTHMM;
        ubounds_CTHMM(i,idx) = ub_CTHMM;
    end

    counter = counter+1;
    disp(['progress: ',num2str(counter),'/',num2str(nRep)])        
end
toc


%% Produce main figure

f = figure('Position',[500 200 600 500]);
tiledlayout(2,2,'TileSpacing','compact','Padding','loose')
% 
x = 1:10:10*nRep;
paramCell = {D1vec,D2vec,k12vec,k21vec};
labelCell = {'D_1 (\mum^2/s)','D_2 (\mum^2/s)','k_{12} (1/s)','k_{21} (1/s)'};

% re-organise data
y_DTHMM = cell2mat(theta_DTHMM);
y_DTHMM = reshape(y_DTHMM,4,[])'; % convert to array
err_DTHMM = cell2mat(error_DTHMM);
err_DTHMM = reshape(err_DTHMM,4,[])'; % convert to array

y_Kinz = cell2mat(theta_Kinz);
y_Kinz = reshape(y_Kinz,4,[])'; % convert to array
err_Kinz = cell2mat(error_Kinz);
err_Kinz = reshape(err_Kinz,4,[])'; % convert to array

y_CTHMM = cell2mat(theta_CTHMM);
y_CTHMM = reshape(y_CTHMM,4,[])'; % convert to array
err_CTHMM = cell2mat(error_CTHMM);
err_CTHMM = reshape(err_CTHMM,4,[])'; % convert to array


for i = 1:4 % do for every parameter

    nexttile
    xconf = [x x(end:-1:1)];

    % Kinz
    % y_iter_Kinz = y_Kinz(:,i)';
    % u_err_iter_Kinz = ubounds_Kinz(:,i)';
    % l_err_iter_Kinz = lbounds_Kinz(:,i)';
    % yconf = [u_err_iter_Kinz l_err_iter_Kinz(end:-1:1)];
    % 
    % p1 = fill(xconf,yconf,'red');
    % p1.FaceColor = [0.9290 0.6940 0.1250];
    % p1.FaceAlpha = 1;
    % p1.EdgeColor = 'none'; 
    % hold on
    % plot(x,y_iter_Kinz,':','Color',[0.9290 0.6940 0.1250]-[0.4 0.4 0.10],'LineWidth',3.0)
    % hold on

    % DTHMM
    y_iter_DTHMM = y_DTHMM(:,i)';
    u_err_iter_DTHMM = ubounds_DTHMM(:,i)';
    l_err_iter_DTHMM = lbounds_DTHMM(:,i)';
    yconf = [u_err_iter_DTHMM l_err_iter_DTHMM(end:-1:1)];
    
    p2 = fill(xconf,yconf,'red');
    p2.FaceColor = [0    0.4470    0.7410];
    p2.FaceAlpha = 0.6;
    p2.EdgeColor = 'none'; 
    hold on
    plot(x,y_iter_DTHMM,'-','Color',[0.3010-0.3    0.7450-0.3    0.9330-0.3],'LineWidth',1.5)
    hold on

    % CTHMM
    y_iter_CTHMM = y_CTHMM(:,i)';
    u_err_iter_CTHMM = ubounds_CTHMM(:,i)';
    l_err_iter_CTHMM = lbounds_CTHMM(:,i)';         
    yconf = [u_err_iter_CTHMM l_err_iter_CTHMM(end:-1:1)];
    p3 = fill(xconf,yconf,'red');
    p3.FaceColor = [1    0.5    0.0];
    p3.FaceAlpha = 0.7;
    p3.EdgeColor = 'none'; 
    hold on
    plot(x,y_iter_CTHMM,'-','Color',[1-0.25 0.5-0.10 0],'LineWidth',1.5)
    hold on

    % ground truth
    p4 = plot(x,paramCell{i},'--k','LineWidth',1.5);

    if i == 1
        legend([p3,p2,p4],{'stHMM','HMM','simulated'},'Location','southwest',...
               'AutoUpdate','off')
        % legend([p3,p1,p2,p4],{'stHMM','BIASD','HMM','simulated'},'Location','northwest',...
        %        'AutoUpdate','off')
    end

    % legend('discrete', 'continuous', 'ground truth')
    ylabel(labelCell{i})
    xlabel('\nu_{1}=\nu_{2}')
    xticks(x)
    xlim([0,x(end)+2])
    if i==3 || i==4
        set(gca,'YScale','log')
    end
    % ylim([0,D1+1])
    % xticks([])
    
    lblCell = cell(1,nRep);
    for j = 1:nRep
    %     lblCell{i} = ['(',num2str(i),',',num2str(i),')'];
        lblCell{j} = [num2str(round(k12vec(j)*tau,2))]; % plot nu_{ij}=k_{ij}*tau
    end
    xticklabels(lblCell)
    box on
end


%% Produce figure of marginalised posterior

f2 = figure('Position',[1000 500 660 350]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')
paramCell = {D1vec,D2vec,k12vec,k21vec};
xlabelCell = {'D_1 (\mum^2/s)','D_2 (\mum^2/s)','k_{12} (1/ms)','k_{21} (1/ms)'};
ylabelCell = {'P(D_1|O)','P(D_2|O)','P(k_{12}|O)','P(k_{21}|O)'};

idxExp = nRep; % choose data set NOTE: this index should be added below to have effect
% [nBins_DTHMM,~,NSoutput_DTHMM,w_DTHMM] = util.NSaxuiliaryPlot(finalSeq_DTHMM,logZ_DTHMM,nLive); % NOTE: extract this info earlier and store it?
% [nBins_CTHMM,~,NSoutput_CTHMM,w_CTHMM] = util.NSaxuiliaryPlot(finalSeq_CTHMM,logZ_CTHMM,nLive); % NOTE: extract this info earlier and store it?

for i = 1:4 % do for every parameter

    nexttile
    % DT-HMM
    [xbar_DTHMM,binCounts_DTHMM] = util.bincounts(nBins_DTHMM, thetaBayes_DTHMM(i), errorBayes_DTHMM(i), NSpos_DTHMM(:,i), w_DTHMM);
    bar(xbar_DTHMM,binCounts_DTHMM,1,'LineStyle','none','FaceColor',[0 0.4470 0.7410],'EdgeColor','none','FaceAlpha',.8) % plot marginalised posterior
    hold on

    % Kinz
    % [xbar_Kinz,binCounts_Kinz] = util.bincounts(nBins_Kinz, thetaBayes_Kinz(i), errorBayes_Kinz(i), NSpos_Kinz(:,i), w_Kinz);
    % bar(xbar_Kinz,binCounts_Kinz,1,'LineStyle','none','FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','none','FaceAlpha',.8) % plot marginalised posterior
    % hold on
    
    % CT-HMM
    [xbar_CTHMM,binCounts_CTHMM] = util.bincounts(nBins_CTHMM, thetaBayes_CTHMM(i), errorBayes_CTHMM(i), NSpos_CTHMM(:,i), w_CTHMM);
    bar(xbar_CTHMM,binCounts_CTHMM,1,'LineStyle','none','FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','none','FaceAlpha',.8) % plot marginalised posterior
    
    % ground truth
    xline(paramCell{i}(idxExp),'--k','LineWidth',1.5)

    % extra labels
    xlabel(xlabelCell{i})
    ylabel(ylabelCell{i})
end

%% add panel names to figure

f1 = openfig('figures/trajectory_main_2025-06-06_4.fig');

ax1 = f1.Children.Children(5);
ax2 = f1.Children.Children(3);
ax3 = f1.Children.Children(2);
% ax4 = f1.Children.Children(4);
ax4 = f1.Children.Children(1);

text(0.02,0.96,{"\bf{a}"},'Parent',ax1,'FontSize',12,'units','normalized');
text(0.02,0.95,{"\bf{b}"},'Parent',ax2,'FontSize',12,'units','normalized');
text(0.02,0.96,{"\bf{c}"},'Parent',ax3,'FontSize',12,'units','normalized');
text(0.02,0.96,{"\bf{d}"},'Parent',ax4,'FontSize',12,'units','normalized');


