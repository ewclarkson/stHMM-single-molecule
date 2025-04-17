function [logZ_error, thetaBayes, errorBayes] = bayesianestimate(NSseq,logZest,nLive,toPlot)

    % 
    % Computes parameter estimates from a nested sampling output sequence
    % 
    % Input:
    % NSSeq = a nested sampling sequence. Structure array with fields
    %         pos (parameter values), postWt (posterior weight) and
    %         logL (log-likelihood)
    % logZest = estimated logarithm of evidence
    % nLive = number of live points in nested sampling run
    % toPlot = 'true' if plotting is desired, 'false' otherwise 
    % 
    % Output:
    % logZ_error = approximate error in estimation of log(evidence)
    % thetaBayes = Bayesian parameter estimates
    % errorBayes = standard deviation of parameter estimates
    % 

    seqLen = length(NSseq); % length of sequence
    
    % Extract data
    data = zeros(seqLen,4); % all coordinate values
    w = zeros(1,seqLen); % corresponding weights
    logLVec = zeros(seqLen,1); % likelihoods
    
    for i = 1:seqLen
        
        data(i,:) = NSseq(i).pos; % parameter values
        w(i) = NSseq(i).postWt; % posterior weights
        logLVec(i) = NSseq(i).logL; % log-likelihood values
    end
    data1 = data(:,1); % 1st coordinate
    data2 = data(:,2); % 2nd coordinate
    data3 = data(:,3); % 3rd coordinate
    data4 = data(:,4); % 4th coordinate
    
    % Compute Bayesian estimates and errors
    logZ_error = sqrt(w*(logLVec-logZest)/nLive);
    thetaBayes = [w*data1,w*data2,w*data3,w*data4]; % posterior means
    M2Bayes = [w*(data1.^2),w*(data2.^2),w*(data3.^2),w*(data4.^2)]; % 2nd moment
    errorBayes = sqrt(M2Bayes-thetaBayes.^2); % standard deviation of estimates

    % Plot marginalised posterior distributions
    if toPlot
        nBins = round(sqrt(seqLen)); % number of bins
    
        figure;
        [xbar,binCounts] = utilB.bincounts(nBins, thetaBayes(1), errorBayes(1), data1, w);
        bar(xbar,binCounts,1) % weighted histogram
        xlabel('D_1')
        ylabel('P(D_1)')
        
        figure;
        [xbar,binCounts] = utilB.bincounts(nBins, thetaBayes(2), errorBayes(2), data2, w);
        bar(xbar,binCounts,1) % weighted histogram
        xlabel('D_2')
        ylabel('P(D_2)')
        
        figure;
        [xbar,binCounts] = utilB.bincounts(nBins, thetaBayes(3), errorBayes(3), data3, w);
        bar(xbar,binCounts,1) % weighted histogram
        xlabel('p_{12}')
        ylabel('P(p_{12})')
        
        figure;
        [xbar,binCounts] = utilB.bincounts(nBins, thetaBayes(4), errorBayes(4), data4, w);
        bar(xbar,binCounts,1) % weighted histogram
        xlabel('p_{21}')
        ylabel('P(p_{21})')
    end
end