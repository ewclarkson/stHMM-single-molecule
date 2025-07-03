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
priorPars1state = {'lognormal',log(0.5*(D1+D2)),log(1.5)};

% nested sampling parameters
nLive = 200; % number of live points
StopRatio = 1E-4; % stop criterion for evidence

% data parameters
nSteps = 100; % trajectory length
nTraj = 10; % number of trajectories

% experimental parameters
sigmaE = 0.010; %0.010; % localisation error std
% M = 100; % number of sub-time intervals for trajectory generation
% mC = 1; % motion blur related. 1: no motion blur
Rmb = 1/6; %1/6; % motion blur coefficient. 1/6: continous shutter mode. 0: no motion blur.
tau = 0.005; % sampling time
exPars = {'tau', tau; 'Rmb', Rmb; 'sigmaE', sigmaE}; % defines experimental parameters that enter models
% NOTE: set Rmb to 1/6 if motion blur is included, otherwise to 0!

% nRows = 5; % number of rows, depending on how many values k12-k21 combinations to test
nStat = 100; % number of analyses of similar datasets

logZcell_DTHMM = cell(1,nRep); % hold estimates of log-evidence
logZstat_DTHMM = zeros(nStat,nRep); % estimates from multiple runs
% logZcell_1state = cell(1,nRep); % hold estimates of log-evidence
logZstat_1state = zeros(nStat,nRep); % hold estimates of log-evidence
logZcell_CTHMM = cell(1,nRep); % hold estimates of log-evidence
logZstat_CTHMM = zeros(nStat,nRep); % estimates from multiple runs

logZcell_error_DTHMM = cell(1,nRep); % hold estimates of log-evidence
% logZcell_error_1state = cell(1,nRep); % hold estimates of log-evidence
logZcell_error_CTHMM = cell(1,nRep); % hold estimates of log-evidence

theta_DTHMM = cell(1,nRep); % to hold errors of parameter estimates
theta_1state = cell(1,nRep); % to hold parameter estimates
theta_CTHMM = cell(1,nRep); % to hold parameter estimates

error_DTHMM = cell(1,nRep); % to hold errors of parameter estimates
error_1state = cell(1,nRep); % to hold errors of parameter estimates
error_CTHMM = cell(1,nRep); % to hold errors of parameter estimates

alpha = 0.025; % follows from chosen confidence level
lbounds_DTHMM = zeros(nRep,4); % lower confidence bounds
ubounds_DTHMM = zeros(nRep,4); % upper confidence bounds
lbounds_1state = zeros(nRep,4); % lower confidence bounds
ubounds_1state = zeros(nRep,4); % upper confidence bounds
lbounds_CTHMM = zeros(nRep,4); % lower confidence bounds
ubounds_CTHMM = zeros(nRep,4); % upper confidence bounds

counter = 0;
for i = 1:nRep
  
    for idy = 1:nStat
        % generate trajectories NOTE: increase the number of sub-time
        % intervals until "convergence"
        data = cell(1,nTraj);
        for k = 1:nTraj
            % data{k} = util.generatetrajectory(nSteps, tau, 1000, D1, D2, k12vec(i), k21vec(i));
            data{k} =  util.generatetrajectory_noisy(nSteps, tau, 1000, 1000, D1, D2, k12vec(i), k21vec(i), sigmaE);
        end
        
        % nested sampling with DT-HMM
        [finalSeq_DTHMM, thetaMLE_DTHMM, logZ_DTHMM] = util.nestedsampling(nLive, StopRatio, priorPars, exPars, data, @util.logl_DTHMM);
        
        % run nested sampling with 1-state diffusion
        [finalSeq_1state,thetaMLE_1state,logZ_1state] = util.nestedsampling_1state(nLive, StopRatio, priorPars1state, exPars, data, @util.logl_1state);
    
        % run nested sampling with CT-HMM
        [finalSeq_CTHMM,thetaMLE_CTHMM,logZ_CTHMM] = util.nestedsampling(nLive, StopRatio, priorPars, exPars, data, @util.logl_CTHMM);
       
        % Compute Bayesian parameter estimates
        [logZ_error_DTHMM, thetaBayes_DTHMM, errorBayes_DTHMM] = util.bayesianestimate(finalSeq_DTHMM,logZ_DTHMM,nLive,0);
        [logZ_error_1state, thetaBayes_1state, errorBayes_1state] = util.bayesianestimate(finalSeq_1state,logZ_1state,nLive,0);
        [logZ_error_CTHMM, thetaBayes_CTHMM, errorBayes_CTHMM] = util.bayesianestimate(finalSeq_CTHMM,logZ_CTHMM,nLive,0);
        
        % store NS output for plotting
        logZcell_DTHMM{i} = logZ_DTHMM;
        logZstat_DTHMM(idy,i) = logZ_DTHMM;
        logZcell_error_DTHMM{i} = logZ_error_DTHMM;
        % logZcell_1state{i} = logZ_1state;
        logZstat_1state(idy,i) = logZ_1state;
        % logZcell_error_DTHMM{i} = logZ_error_Kinz;
        logZcell_CTHMM{i} = logZ_CTHMM;
        logZstat_CTHMM(idy,i) = logZ_CTHMM;
        logZcell_error_CTHMM{i} = logZ_error_CTHMM;
    
        theta_DTHMM{i} = thetaBayes_DTHMM;
        error_DTHMM{i} = errorBayes_DTHMM;
        % theta_1state{i} = thetaBayes_1state;
        % error_1state{i} = errorBayes_1state;
        theta_CTHMM{i} = thetaBayes_CTHMM;
        error_CTHMM{i} = errorBayes_CTHMM;
        
        [nBins_DTHMM,~,NSpos_DTHMM,w_DTHMM] = util.NSaxuiliaryPlot(finalSeq_DTHMM,logZ_DTHMM,nLive);
        % [nBins_1state,~,NSpos_1state,w_1state] = util.NSaxuiliaryPlot(finalSeq_1state,logZ_1state,nLive);
        [nBins_CTHMM,~,NSpos_CTHMM,w_CTHMM] = util.NSaxuiliaryPlot(finalSeq_CTHMM,logZ_CTHMM,nLive);
        
        for idx = 1:4 % NOTE: put this loop inside confidenceinterval.m
            [xbar_DTHMM,binCounts_DTHMM] = util.bincounts(nBins_DTHMM, thetaBayes_DTHMM(idx), errorBayes_DTHMM(idx), NSpos_DTHMM(:,idx), w_DTHMM);
            % [xbar_Kinz,binCounts_1state] = util.bincounts(nBins_1state, thetaBayes_1state(idx), errorBayes_1state(idx), NSpos_1state(:,idx), w_1state);
            [xbar_CTHMM,binCounts_CTHMM] = util.bincounts(nBins_CTHMM, thetaBayes_CTHMM(idx), errorBayes_CTHMM(idx), NSpos_CTHMM(:,idx), w_CTHMM);
            [lb_DTHMM,ub_DTHMM] = util.ci(xbar_DTHMM,binCounts_DTHMM,alpha,alpha);
            % [lb_1state,ub_1state] = util.ci(xbar_1state,binCounts_1state,alpha,alpha);
            [lb_CTHMM,ub_CTHMM] = util.ci(xbar_CTHMM,binCounts_CTHMM,alpha,alpha);
            lbounds_DTHMM(i,idx) = lb_DTHMM;
            ubounds_DTHMM(i,idx) = ub_DTHMM;
            % lbounds_1state(i,idx) = lb_1state;
            % ubounds_1state(i,idx) = ub_1state;
            lbounds_CTHMM(i,idx) = lb_CTHMM;
            ubounds_CTHMM(i,idx) = ub_CTHMM;
        end
    end
    counter = counter+1;
    disp(['progress: ',num2str(counter),'/',num2str(nRep)])        
end
toc

%% Plot model probabilities

f = figure('Position',[500 300 600 300]);
% tiledlayout(2,2,'TileSpacing','loose','Padding','loose')
% 
x = 1:10:10*nRep;

% Re-organise data
% logZarr_DTHMM = cell2mat(logZcell_DTHMM);
% logZarr_CTHMM = cell2mat(logZcell_CTHMM);
% logZarr_error_DTHMM = cell2mat(logZcell_error_DTHMM);
% logZarr_error_CTHMM = cell2mat(logZcell_error_DTHMM);
logZarr_DTHMM = mean(logZstat_DTHMM); % mean value of runs
logZarr_1state = mean(logZstat_1state); % mean value of runs
logZarr_CTHMM = mean(logZstat_CTHMM); % mean value of runs

stHMMprob = zeros(1,7);
onestateprob = zeros(1,7);
HMMprob = zeros(1,7);
% mProb_upper = zeros(1,7);
% mProb_lower = zeros(1,7);
for i = 1:nRep
    stHMMprob(i) = exp(logZarr_CTHMM(i)-util.logsumexp2(util.logsumexp2(logZarr_CTHMM(i),logZarr_DTHMM(i)),logZarr_1state(i)));
    onestateprob(i) = exp(logZarr_1state(i)-util.logsumexp2(util.logsumexp2(logZarr_CTHMM(i),logZarr_DTHMM(i)),logZarr_1state(i)));
    HMMprob(i) = exp(logZarr_DTHMM(i)-util.logsumexp2(util.logsumexp2(logZarr_CTHMM(i),logZarr_DTHMM(i)),logZarr_1state(i))); % follows from normalisation

    % mProb_upper(i) = exp(logZarr_CTHMM(i)+logZarr_error_CTHMM(i)-util.logsumexp2...
    %     (logZarr_CTHMM(i)+logZarr_error_CTHMM(i),logZarr_DTHMM(i)-logZarr_error_DTHMM(i)));
    % mProb_lower(i) = exp(logZarr_CTHMM(i)-logZarr_error_CTHMM(i)-util.logsumexp2...
    %     (logZarr_CTHMM(i)-logZarr_error_CTHMM(i),logZarr_DTHMM(i)+logZarr_error_DTHMM(i)));
end

% xconf = [x x(end:-1:1)];
% yconf = [mProb_upper mProb_lower(end:-1:1)];
% 
% p = fill(xconf,yconf,'red');
% p.FaceColor = [0.3010 0.7450 0.9330];
% p.FaceAlpha = 0.8;
% p.EdgeColor = 'none'; 
% hold on
plot(x,stHMMprob,'Color',[0 0.4470 0.7410],'LineWidth',1.5)
hold on
plot(x,HMMprob,'Color',[0.8500    0.3250    0.0980],'LineWidth',1.5)
hold on
plot(x,onestateprob,'Color',[0.9290    0.6940    0.1250],'LineWidth',1.5)
legend('stHMM','HMM','one-state', 'Location','southwest')

% legend('discrete', 'continuous', 'ground truth')
ylabel('model probability')
% ylim([0.5,1])
xlabel('\nu_1=\nu_2')
xticks(x)
xlim([0,x(end)+2])

lblCell = cell(1,nRep);
for j = 1:nRep
%     lblCell{i} = ['(',num2str(i),',',num2str(i),')'];
    lblCell{j} = [num2str(round(k12vec(j)*tau,2))];
end
xticklabels(lblCell)
box on



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



%% plot the model selection results

f1 = openfig('figures/model_probabilities_D1=10_2025-07-01.fig');
ax1 = gca;

f2 = openfig('figures/model_probabilities_D1=3_2025-06-24.fig');
% f1.Position = [500 200 600 500];
ax2 = gca;

f3 = figure('Position',[500 100 600 350]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

ax1c = copyobj(ax1,tl);
ax1c.Layout.Tile = 1;
ax1c.Layout.TileSpan = [1 1];

ax2c = copyobj(ax2,tl);
ax2c.Layout.Tile = 2;
ax2c.Layout.TileSpan = [1 1];

legend('stHMM','HMM','one-state')