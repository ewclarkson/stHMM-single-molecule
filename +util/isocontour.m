function [xbar1,xbar2,binCounts] = isocontour(NSseq, nBins, thetaBayes, stdTheta, pair)

    % 
    % Plot iso-contours of a sampled distribution
    % 
    % Input:
    % NSseq = sequence of posterior weights and parameter positions sampled
    %         via the nested sampling algorithm
    % nBins = number of bins to use for the weighted histogram
    % thetaBayes = mean values (estimates) of parameters of interest
    % stdTheta = standard deviations of estimates
    % pair = two parameters of interest, on the form [1,3]
    % 
    % Output:
    % xbar = vector of bin centra
    % binCounts = vector of bin counts in each bin
    % 
    
    seqLen = length(NSseq); % length of sequence

    % extract data
    data = zeros(seqLen,4); % all coordinate values
    w = zeros(1,seqLen); % corresponding weights
    logLVec = zeros(seqLen,1); % likelihoods

    for i = 1:seqLen

        data(i,:) = NSseq(i).pos; % parameter values
        w(i) = NSseq(i).postWt; % posterior weights
        logLVec(i) = NSseq(i).logL; % log-likelihood values
    end
    data1 = data(:,pair(1)); % 1st coordinate
    data2 = data(:,pair(2)); % 2nd coordinate

    % for the first parameter
    edges1 = zeros(1,nBins+1); % bin edges
    binStart1 = thetaBayes(1)-5*stdTheta(1);
    binStop1 = thetaBayes(1)+5*stdTheta(1);
    binWidth1 = (binStop1-binStart1)/nBins;

    % for the second parameter
    edges2 = zeros(1,nBins+1); % bin edges
    binStart2 = thetaBayes(2)-5*stdTheta(2);
    binStop2 = thetaBayes(2)+5*stdTheta(2);
    binWidth2 = (binStop2-binStart2)/nBins;
    
    % add bin edges
    for k = 1:nBins+1
        edges1(k) = binStart1 + (k-1)*binWidth1;
        edges2(k) = binStart2 + (k-1)*binWidth2;
    end
    
    % fill bins and add weights
    binCounts = zeros(nBins,nBins);
    for i = 1:seqLen % go through all points
        for j = 1:nBins % find the bin to contain this point
            for k = 1:nBins
           
                if (edges1(j) < data1(i)) && (data1(i) <= edges1(j+1)) && (edges2(k) < data2(i)) && (data2(i) <= edges2(k+1))
        
                    binCounts(j,k) = binCounts(j,k) + w(i);
                end
            end
        end
    end

    binCounts = binCounts/(binWidth1*binWidth2); % normalisation
    % trapz(binWidth1,trapz(binWidth2,binCounts,2)); % test normalisation
    binStop1 = edges1(end);
    binStop2 = edges2(end);
    xbar1 = binStart1+binWidth1:binWidth1:binStop1; % bin centra
    xbar2 = binStart2+binWidth2:binWidth2:binStop2; % bin centra
end



