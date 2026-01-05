function [deadPoints,thetaMLE,logZ] = ...
    nestedsampling_1state(nLive, stopRatio, priorPars, exPars, data, loglfun)

    % 
    % Perform the nested sampling algorithm in the special case of only 1
    % parameter.
    % 
    % Input:
    % nLive     - number of live points to be used
    % StopRatio - stopping criterion as the fractional change in estimated 
    %             evidence falls below this value
    % priorPars - array of upper and lower boundaries for uniform priors
    % exPars    - cell array of experimental parameters.
    %             exPars = {'tau', 5; 'Rmb', 1/6; 'sigmaE', 0.1}
    %             where 'tau' is the sampling time, 'Rmb' is the motion blur
    %             coefficient and 'sigmaE' is the localisation error std.
    %             NOTE: simply set Rmb=sigmaE=0 in case no nosie corrections
    %             are to be used.
    % data      - cell array containing trajectories of displacements
    % loglfun   - function handle to chosen likelihood function
    % 
    % Output:
    % deadPoints - samled coordinates with posterior weights
    % thetaMLE   - maximum likelihood estimate of parameters
    % logZ       - estimated logarithm of evidence
    % 
    % Dependencies:
    % logsumexp2.m
    % drawlivepoint.m
    % 

    % ----------------- set up live points ---------------------------

    % Initialise variables
    logStopRatio = log(stopRatio); % logarithm of evidence ratio
    logZ = log(0); % initialise evidence
    logZratio = log(Inf); % initialise log-ratio of estimated evidence
    j = 1; % iteration index
 
    % Sample live points uniformly in unit hypercube
    uMat = rand(nLive,1); % (point index, parameter index) 
    
    % Compute positions in real parameter space
    posMat = zeros(nLive,1);
    for k = 1:1
        if strcmp(priorPars{k,1},'uniform')
            posMat(:,k) = unifinv(uMat(:,k),priorPars{k,2},priorPars{k,3});
        
        elseif strcmp(priorPars(k,1),'lognormal')
            posMat(:,k) = logninv(uMat(:,k),priorPars{k,2},priorPars{k,3});
        end
        % posMat(:,k) = unifinv(uMat(:,k),priorPars(k,1),priorPars(k,2));
    end

    % Assign parameter values to each live point
    for i = 1:nLive

       u = uMat(i,:); 
       y = posMat(i,:);    
    
       livePoints(i).upos = u; % assign values to points
       livePoints(i).logL = loglfun(y, exPars, data); % compute and store log-likelihood
       livePoints(i).pos = y; % assign values to points
    end
    
    % ----------------------- begin algorithm ----------------------

    while logZratio > logStopRatio % do until the stop criterion is met, i.e. until Z has been calculated
    
        % Identify the worst point
        [logLworst, ind_worst]=min([livePoints.logL]);
        
        % Compute evidence Z up to the jth iteration
        if j==1 % do for the first iteration
            logWeight = -log(nLive+1);
        else % do for all other iterations
            logWeight = logWeight+log(nLive)-log(nLive+1);
        end
        logZ = util.logsumexp2(logZ, logLworst+logWeight);

        % --------------- debugging -----------------------
        % if logninv(livePoints(ind_worst).upos(1),priorPars{1,2},priorPars{1,3}) ~= livePoints(ind_worst).pos(1) % debugging
        % 
        %     keyboard
        % end
        % --------------- end of debugging -----------------
        
        % Sample a new point 
        [uposNew, logLnew, posNew] = util.drawlivepoint_1state(livePoints, logLworst, exPars, data, priorPars, loglfun); % new point
        
        % Store dead point
        deadPoints(j).pos = livePoints(ind_worst).pos;
        deadPoints(j).logL = logLworst;
        deadPoints(j).logWeight = logWeight;

        % Replace dead point
        livePoints(ind_worst).upos = uposNew; % set position of new point
        livePoints(ind_worst).logL = logLnew; % set log-likelihood of new point
        livePoints(ind_worst).pos = posNew; % set position of new point
        j = j+1; % update iteration number
        
        % Compute remaining evidence of points, at iteration j
        logLremain = livePoints(1).logL;
        for i = 2:nLive % do for every live point
            logLremain = util.logsumexp2(logLremain,livePoints(i).logL);
        end
        logZremain = logWeight+logLremain; % remaining log-evidence
    
        logZratio = logZremain-logZ;
    end
    
    % Begin final correction
    logZ = util.logsumexp2(logZ, logZremain);

    % Compute all posterior weights    
    for i = 1:length(deadPoints) % do for every dead point
        
        deadPoints(i).postWt = exp(deadPoints(i).logWeight+...
            deadPoints(i).logL-logZ);     
    end
    
    for i = 1:nLive % do for every live point left
        
        deadPoints(end+1).postWt = exp(logWeight+...
            livePoints(i).logL-logZ);

        deadPoints(end).pos = livePoints(i).pos;
        deadPoints(end).logWeight = logWeight;
        deadPoints(end).logL = livePoints(i).logL;
    end

    % Find MLE of model parameters theta
    [~, ind_max]=max([livePoints.logL]);
    thetaMLE = livePoints(ind_max).upos; % in transformed space
end


