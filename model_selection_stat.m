%% Do model selection for stHMM, HMM and one-state diffusion models

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

% set priors (same for all nRep datasets)
priorPars = {'lognormal',log(D1),log(1.5); 'lognormal',log(D2),log(1.5); 'lognormal',log(200),log(3.5); 'lognormal',log(200),log(3.5)}; 
priorPars1state = {'lognormal',log(0.5*(D1+D2)),log(1.5)}; % 1-state prior for diffusion constant

% nested sampling parameters
nLive = 200; % number of live points
StopRatio = 1E-4; % stop criterion for evidence

% data parameters
nSteps = 100; % trajectory length
nTraj = 10; % number of trajectories

% experimental parameters
sigmaE = 0.010; %0.010; % localisation error std
Rmb = 1/6; %1/6; % motion blur coefficient. 1/6: continous shutter mode. 0: no motion blur.
tau = 0.005; % sampling time
exPars = {'tau', tau; 'Rmb', Rmb; 'sigmaE', sigmaE}; % defines experimental parameters that enter models
% NOTE: set Rmb to 1/6 if motion blur is included, otherwise set it to 0
nStat = 100; % number of analyses of similar datasets

logZcell_DTHMM = cell(1,nRep); % hold estimates of log-evidence
logZstat_DTHMM = zeros(nStat,nRep); % estimates from multiple runs
logZstat_1state = zeros(nStat,nRep); % hold estimates of log-evidence
logZcell_CTHMM = cell(1,nRep); % hold estimates of log-evidence
logZstat_CTHMM = zeros(nStat,nRep); % estimates from multiple runs

logZcell_error_DTHMM = cell(1,nRep); % hold estimates of log-evidence
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
        % generate trajectories
        data = cell(1,nTraj);
        for k = 1:nTraj
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
        logZstat_1state(idy,i) = logZ_1state;
        logZcell_CTHMM{i} = logZ_CTHMM;
        logZstat_CTHMM(idy,i) = logZ_CTHMM;
        logZcell_error_CTHMM{i} = logZ_error_CTHMM;
    
        theta_DTHMM{i} = thetaBayes_DTHMM;
        error_DTHMM{i} = errorBayes_DTHMM;
        theta_CTHMM{i} = thetaBayes_CTHMM;
        error_CTHMM{i} = errorBayes_CTHMM;
        
        [nBins_DTHMM,~,NSpos_DTHMM,w_DTHMM] = util.NSaxuiliaryPlot(finalSeq_DTHMM,logZ_DTHMM,nLive);
        [nBins_CTHMM,~,NSpos_CTHMM,w_CTHMM] = util.NSaxuiliaryPlot(finalSeq_CTHMM,logZ_CTHMM,nLive);
        
        for idx = 1:4
            [xbar_DTHMM,binCounts_DTHMM] = util.bincounts(nBins_DTHMM, thetaBayes_DTHMM(idx), errorBayes_DTHMM(idx), NSpos_DTHMM(:,idx), w_DTHMM);
            [xbar_CTHMM,binCounts_CTHMM] = util.bincounts(nBins_CTHMM, thetaBayes_CTHMM(idx), errorBayes_CTHMM(idx), NSpos_CTHMM(:,idx), w_CTHMM);
            [lb_DTHMM,ub_DTHMM] = util.ci(xbar_DTHMM,binCounts_DTHMM,alpha,alpha);
            [lb_CTHMM,ub_CTHMM] = util.ci(xbar_CTHMM,binCounts_CTHMM,alpha,alpha);
            lbounds_DTHMM(i,idx) = lb_DTHMM;
            ubounds_DTHMM(i,idx) = ub_DTHMM;
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
x = 1:10:10*nRep;

% Re-organise data
logZarr_DTHMM = mean(logZstat_DTHMM); % mean value of runs
logZarr_1state = mean(logZstat_1state); % mean value of runs
logZarr_CTHMM = mean(logZstat_CTHMM); % mean value of runs

stHMMprob = zeros(1,7);
onestateprob = zeros(1,7);
HMMprob = zeros(1,7);
for i = 1:nRep
    stHMMprob(i) = exp(logZarr_CTHMM(i)-util.logsumexp2(util.logsumexp2(logZarr_CTHMM(i),logZarr_DTHMM(i)),logZarr_1state(i)));
    onestateprob(i) = exp(logZarr_1state(i)-util.logsumexp2(util.logsumexp2(logZarr_CTHMM(i),logZarr_DTHMM(i)),logZarr_1state(i)));
    HMMprob(i) = exp(logZarr_DTHMM(i)-util.logsumexp2(util.logsumexp2(logZarr_CTHMM(i),logZarr_DTHMM(i)),logZarr_1state(i))); % follows also from normalisation
end

plot(x,stHMMprob,'Color',[0 0.4470 0.7410],'LineWidth',1.5)
hold on
plot(x,HMMprob,'Color',[0.8500    0.3250    0.0980],'LineWidth',1.5)
hold on
plot(x,onestateprob,'Color',[0.9290    0.6940    0.1250],'LineWidth',1.5)
legend('stHMM','HMM','one-state', 'Location','southwest')

ylabel('model probability')
% ylim([0.5,1])
xlabel('\nu_1=\nu_2')
xticks(x)
xlim([0,x(end)+2])

lblCell = cell(1,nRep);
for j = 1:nRep
    lblCell{j} = [num2str(round(k12vec(j)*tau,2))];
end
xticklabels(lblCell)
box on


%% plot the model selection results. NOTE: first save sub-figures above.

f1 = openfig('figures/model_probabilities_D1=10_2025-07-01.fig');
ax1 = gca;

f2 = openfig('figures/model_probabilities_D1=3_2025-06-24.fig');
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
