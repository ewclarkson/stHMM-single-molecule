%% Compare continous-time HMM with and without noise sources for parameter inference

% assume D1 > D2
tic
nRep = 4; % try this many values of k12,k21

% model parameters
D1 = 10; % diffusion constant of state 1
D2 = 1; % diffusion constant of state 2
D1vec = linspace(D1,D1,nRep);
D2vec = linspace(D2,D2,nRep);
k12_min_gen = 1E-2/5;
k12_max_gen = 1/5;
k21_min_gen = 1E-2/5;
k21_max_gen = 1/5;
% k12vec = logspace(log10(k12_min_gen),log10(k12_max_gen),nRep);
% k21vec = logspace(log10(k21_min_gen),log10(k21_max_gen),nRep);
k12 = 0.02; % NOTE: what value do we want here?
k21 = k12;
k12vec = linspace(k12,k12,nRep);
k21vec = linspace(k21,k21,nRep);

% nested sampling parameters
nLive = 100; % number of live points
StopRatio = 1E-4; % stop criterion for evidence

% data parameters
nSteps = 300; % trajectory length
nTraj = 10; % number of trajectories

% experimental parameters
% sigmaE = 0; % localisation error std
% M = 100; % number of sub-time intervals for trajectory generation
% mC = 1; % motion blur related. 1: no motion blur
% Rmb = 0; % motion blur coefficient. 1/6: continous shutter mode. 0: no motion blur.
tau = 5; % sampling time
sigmaE_min = 0; % NOTE: is this formula correct?
sigmaE_max = 0.2*sqrt(2*D2*tau);
sigmaEvec = linspace(sigmaE_min,sigmaE_max,nRep);

% nRows = 5; % number of rows, depending on how many values k12-k21 combinations to test

logZ_DTHMM_all = zeros(1,nRep); % hold estimates of log-evidence
logZ_CTHMM_all = zeros(1,nRep); % hold estimates of log-evidence
logZ_error_DTHMM_all = zeros(1,nRep);
logZ_error_CTHMM_all = zeros(1,nRep);

theta_DTHMM = cell(1,nRep); % to hold errors of parameter estimates
theta_CTHMM = cell(1,nRep); % to hold parameter estimates
error_DTHMM = cell(1,nRep); % to hold errors of parameter estimates
error_CTHMM = cell(1,nRep); % to hold errors of parameter estimates

alpha = 0.025; % follows from chosen confidence level
lbounds_DTHMM = zeros(nRep,4); % lower confidence bounds
ubounds_DTHMM = zeros(nRep,4); % upper confidence bounds
% lbounds_Kinz = zeros(nRep,4); % lower confidence bounds
% ubounds_Kinz = zeros(nRep,4); % upper confidence bounds
lbounds_CTHMM = zeros(nRep,4); % lower confidence bounds
ubounds_CTHMM = zeros(nRep,4); % upper confidence bounds

counter = 0;
for i = 1:nRep
  
    % generate trajectories NOTE: increase the number of sub-time
    % intervals until "convergence"

    exPars = {'tau', tau; 'Rmb', 1/6; 'sigmaE', sigmaEvec(i)}; % defines experimental parameters that enter models

    data = cell(1,nTraj);
    for k = 1:nTraj
        % data{k} = util.generatetrajectory(k12,k21,D1,D2,tau,nSteps); % no noise
        data{k} = util.noisyBrownian2state2D(nSteps, tau, 100, 100, D1, D2, k12, k21, sigmaEvec(i)); % with noise
    end
     
    % set priors
    D1_max = 20; % known upper limit of D1
    D1_min = 2; % known lower limit of D1
    D2_max = 5; % known upper limit of D2
    D2_min = 0.01; % known lower limit of D2
    k12_min = 1E-3;
    k12_max = 3*k12_max_gen;
    k21_min = 1E-3;
    k21_max = 3*k21_max_gen;
    % priorLimits = [D1_min,D1_max;D2_min,D2_max;k12_min,k12_max;k21_min,k21_max];
    priorPars = {'uniform',D1_min,D1_max; 'uniform',D2_min,D2_max; 'uniform',k12_min,k12_max; 'uniform',k21_min,k21_max}; 
    
    % nested sampling with DT-HMM % NOTE: update the name of this
    [finalSeq_DTHMM, thetaMLE_DTHMM, logZ_DTHMM] = util.nestedsampling(nLive, StopRatio, priorPars, exPars, data, @util.logl_CTHMM_noisecorr);
    
    % run nested sampling with CT-HMM
    [finalSeq_CTHMM,thetaMLE_CTHMM,logZ_CTHMM] = util.nestedsampling(nLive, StopRatio, priorPars, exPars, data, @util.logl_CTHMM);
   
    % Compute Bayesian parameter estimates
    [logZ_error_DTHMM, thetaBayes_DTHMM, errorBayes_DTHMM] = util.bayesianestimate(finalSeq_DTHMM,logZ_DTHMM,nLive,0);
    [logZ_error_CTHMM, thetaBayes_CTHMM, errorBayes_CTHMM] = util.bayesianestimate(finalSeq_CTHMM,logZ_CTHMM,nLive,0);
    
    % store NS output for plotting
    theta_DTHMM{i} = thetaBayes_DTHMM;
    error_DTHMM{i} = errorBayes_DTHMM;
    theta_CTHMM{i} = thetaBayes_CTHMM;
    error_CTHMM{i} = errorBayes_CTHMM;

    logZ_DTHMM_all(i) = logZ_DTHMM;
    logZ_CTHMM_all(i) = logZ_CTHMM;
    
    [nBins_DTHMM,~,NSpos_DTHMM,w_DTHMM] = util.NSaxuiliaryPlot(finalSeq_DTHMM,logZ_DTHMM,nLive);
    [nBins_CTHMM,~,NSpos_CTHMM,w_CTHMM] = util.NSaxuiliaryPlot(finalSeq_CTHMM,logZ_CTHMM,nLive);
    
    % for idx = 1:4 % NOTE: put this loop inside confidenceinterval.m
    %     [xbar_DTHMM,binCounts_DTHMM] = util.bincounts(nBins_DTHMM, thetaBayes_DTHMM(idx), errorBayes_DTHMM(idx), NSpos_DTHMM(:,idx), w_DTHMM);
    %     % [xbar_Kinz,binCounts_Kinz] = util.bincounts(nBins_Kinz, thetaBayes_Kinz(idx), errorBayes_Kinz(idx), NSpos_Kinz(:,idx), w_Kinz);
    %     [xbar_CTHMM,binCounts_CTHMM] = util.bincounts(nBins_CTHMM, thetaBayes_CTHMM(idx), errorBayes_CTHMM(idx), NSpos_CTHMM(:,idx), w_CTHMM);
    %     [lb_DTHMM,ub_DTHMM] = util.ci(xbar_DTHMM,binCounts_DTHMM,alpha,alpha);
    %     % [lb_Kinz,ub_Kinz] = util.ci(xbar_Kinz,binCounts_Kinz,alpha,alpha);
    %     [lb_CTHMM,ub_CTHMM] = util.ci(xbar_CTHMM,binCounts_CTHMM,alpha,alpha);
    %     lbounds_DTHMM(i,idx) = lb_DTHMM;
    %     ubounds_DTHMM(i,idx) = ub_DTHMM;
    %     % lbounds_Kinz(i,idx) = lb_Kinz;
    %     % ubounds_Kinz(i,idx) = ub_Kinz;
    %     lbounds_CTHMM(i,idx) = lb_CTHMM;
    %     ubounds_CTHMM(i,idx) = ub_CTHMM;
    % end

    counter = counter+1;
    disp(['progress: ',num2str(counter),'/',num2str(nRep)])        
end
toc


%% Produce main figure

f = figure('Position',[1000 500 290 300]);
% tiledlayout(2,2,'TileSpacing','loose','Padding','loose')
% 
x = 1:10:10*nRep;

% Re-organise data
logBF_all = logZ_DTHMM_all-logZ_CTHMM_all;
logBF_upper = logBF_all+logZ_error_DTHMM+logZ_error_CTHMM;
logBF_lower = logBF_all-logZ_error_DTHMM-logZ_error_CTHMM;
% paramCell = {D1vec,D2vec,k12vec,k21vec};
% labelCell = {'D_1 (\mum^2/s)','D_2 (\mum^2/s)','k_{12} (1/ms)','k_{21} (1/ms)'};

% for i = 1:4 % do for every parameter
% 
%     nexttile
    xconf = [x x(end:-1:1)];

    % log(BF)
    % yconf = [logZ_DTHMM_all+logZ_error_DTHMM_all logZ_DTHMM_all(end:-1:1)-logZ_error_DTHMM_all(end:-1:1)];
    yconf = [logBF_upper logBF_lower(end:-1:1)];
    
    p = fill(xconf,yconf,'red');
    p.FaceColor = [0.3010 0.7450 0.9330];
    p.FaceAlpha = 0.8;
    p.EdgeColor = 'none'; 
    hold on
    plot(x,logBF_all,'Color',[0 0.4470 0.7410],'LineWidth',1.5)

    % legend('discrete', 'continuous', 'ground truth')
    ylabel('log(BF)')
    xlabel('\sigma_{E}')
    xticks(x)
    xlim([0,x(end)+2])
    
    lblCell = cell(1,nRep);
    for j = 1:nRep
    %     lblCell{i} = ['(',num2str(i),',',num2str(i),')'];
        lblCell{j} = [num2str(round(sigmaEvec(j),2))];
    end
    xticklabels(lblCell)
    box on
% end


%% Produce figure of marginalised posterior

f2 = figure('Position',[1000 500 580 350]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')
paramCell = {D1vec,D2vec,k12vec,k21vec};
xlabelCell = {'D_1 (\mum^2/s)','D_2 (\mum^2/s)','k_{12} (1/ms)','k_{21} (1/ms)'};
ylabelCell = {'P(D_1|O)','P(D_2|O)','P(k_{12}|O)','P(k_{21}|O)'};

idxExp = nRep; % choose data set NOTE: this index should be added below to have effect
% [nBins_DTHMM,~,NSoutput_DTHMM,w_DTHMM] = util.NSaxuiliaryPlot(finalSeq_DTHMM,logZ_DTHMM,nLive); % NOTE: extract this info earlier and store it?
% [nBins_CTHMM,~,NSoutput_CTHMM,w_CTHMM] = util.NSaxuiliaryPlot(finalSeq_CTHMM,logZ_CTHMM,nLive); % NOTE: extract this info earlier and store it?

for i = 1:4 % do for every parameter

    nexttile
    % CT-HMM
    [xbar_CTHMM,binCounts_CTHMM] = util.bincounts(nBins_CTHMM, thetaBayes_CTHMM(i), errorBayes_CTHMM(i), NSpos_CTHMM(:,i), w_CTHMM);
    bar(xbar_CTHMM,binCounts_CTHMM,1,'LineStyle','none','FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','none','FaceAlpha',.9) % plot marginalised posterior
    hold on

    % DT-HMM
    [xbar_DTHMM,binCounts_DTHMM] = util.bincounts(nBins_DTHMM, thetaBayes_DTHMM(i), errorBayes_DTHMM(i), NSpos_DTHMM(:,i), w_DTHMM);
    bar(xbar_DTHMM,binCounts_DTHMM,1,'LineStyle','none','FaceColor',[0.4660 0.6740 0.1880],'EdgeColor','none','FaceAlpha',.9) % plot marginalised posterior
    
    % ground truth
    xline(paramCell{i}(idxExp),'--k','LineWidth',1.5)

    % extra labels
    xlabel(xlabelCell{i})
    ylabel(ylabelCell{i})
end





