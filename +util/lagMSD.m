
function [MSDlagXY,MSDlagX,MSDlagY] = lagMSD(data,maxLen)
    
    % Input:
    % data = cell array of doubles
    % maxLen = number of positions in the longest trajectory

    numTraj = length(data);
    dispMatXY = zeros(maxLen-1,numTraj);
    dispMatX = zeros(maxLen-1,numTraj);
    dispMatY = zeros(maxLen-1,numTraj);

    for idv = 1:numTraj % trajectory index
        traj = data{idv}; 
        len = length(traj);
        % lagMSDvec = zeros(1,len-1); 
        for idx = 1:len-1 % lag index
            lagMSDxy = 0; % initialise
            lagMSDx = 0;
            lagMSDy = 0;
            for idy = 1:len-idx % position index
                lagMSDxy = lagMSDxy + sum((traj(idy+idx,:)-traj(idy,:)).^2);
                lagMSDx = lagMSDx + (traj(idy+idx,1)-traj(idy,1))^2;
                lagMSDy = lagMSDy + (traj(idy+idx,2)-traj(idy,2))^2;
            end
            dispMatXY(idx,idv) = 1/(len-idx)*lagMSDxy;
            dispMatX(idx,idv) = 1/(len-idx)*lagMSDx;
            dispMatY(idx,idv) = 1/(len-idx)*lagMSDy;
        end
    end
    
    MSDlagXY = zeros(maxLen-1,1);
    MSDlagX = zeros(maxLen-1,1);
    MSDlagY = zeros(maxLen-1,1);
    for idx = 1:maxLen-1
        MSDlagXY(idx) = sum(dispMatXY(idx,:))/nnz(dispMatXY(idx,:));
        MSDlagX(idx) = sum(dispMatX(idx,:))/nnz(dispMatX(idx,:));
        MSDlagY(idx) = sum(dispMatY(idx,:))/nnz(dispMatY(idx,:));
    end
end




