function [lBound,rBound] = ci(binCentra,binCounts,alpha,beta)
    % 
    % Compute confidence levels from given histogram
    % 
    % Input:
    % binCentra - centra of histogram bins
    % binCounts - normalised bin counts of histogram
    % alpha - lower confidence level, e.g. 0.05
    % beta - upper confidence level, usually set = alpha
    % 
    % Output:
    % lBound - left bound, i.e. x-value (bin centre) corresponding to alpha
    % rBound - right bound, i.e. x-value (bin centre) corresponding to beta
    % 

    % use the example below to test the code
    % h = 0.1
    % binCounts = [0.4,0.7,1.3,1.5,1.9,1.3,1.1,0.9,0.6,0.3];
    % binCentra = 0.2:h:1.1;
    
    % sum(h*binCounts) % test (approximate) normalisation
    
    % compute bin spacing
    h = binCentra(2)-binCentra(1); % assumes uniform spacing

    % find the lower confidence level
    CDFsum = cumsum(h*binCounts);
    % alpha = 0.1;
    indS = find(CDFsum<=alpha,1,'last');
    lBound = binCentra(indS); % lower confidence level
    
    % find the upper confidence level
    CDFsum2 = cumsum(h*flip(binCounts));
    % beta = alpha;
    indS2 = find(CDFsum2<=beta,1,'last');
    rBound = binCentra(end-indS2); % upper confidence level

end



