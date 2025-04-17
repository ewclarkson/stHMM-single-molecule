function [xbar,binCounts] = bincounts(nBins, thetaBayes, stdTheta, data, w)

    % 
    % Bin data for making a weighted histogram of marginalised
    % probability distributions
    % 
    % Input:
    % nBins = number of bins
    % thetaBayes = weighted mean of data
    % stdTheta = standard deviation of estimate
    % data = vector of data to be binned
    % 
    % Output:
    % xbar = vector of bin centra
    % binCounts = vector of bin counts in each bin
    % 

    % find a reasonable x-interval
    seqLen = length(data);

    edges = zeros(1,nBins+1); % bin edges
    binStart = thetaBayes(1)-5*stdTheta(1);
    binStop = thetaBayes(1)+5*stdTheta(1);
    binWidth = (binStop-binStart)/nBins;
    
    for k = 1:nBins+1
        edges(k) = binStart + (k-1)*binWidth;
    end
    
    % fill bins and add weights
    binCounts = zeros(1,nBins);
    for i = 1:seqLen % go through all points
        for j = 1:nBins % find the bin to contain this point
    
            if (edges(j) < data(i)) && (data(i) <= edges(j+1))
    
                binCounts(j) = binCounts(j) + w(i);
            end
        end
    end
    binCounts = binCounts/binWidth; % normalisation
    binStop = edges(end);
    xbar = binStart+binWidth:binWidth:binStop; % bin centra
end

