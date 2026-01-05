function [nBins,thetaBayes,data,w] = NSaxuiliaryPlot(NSseq,logZest,nLive)

    seqLen = length(NSseq); % length of sequence
    nBins = round(sqrt(seqLen)); % number of bins for later plotting
    
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

end
