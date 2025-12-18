function [MSD_xy,MSD_x,MSD_y] = MSD(data,maxLen)
    
    % Input:
    % data = cell array of doubles
    % maxLen = number of positions in the longest trajectory

    numTraj = length(data);
    dispMatXY = zeros(maxLen-1,numTraj);
    dispMatX = zeros(maxLen-1,numTraj);
    dispMatY = zeros(maxLen-1,numTraj);
    
    for idx = 1:numTraj
    
        traj = data{idx};
        len = length(traj);
        disp = traj(2:end,:)-repmat(traj(1,:),len-1,1); % subtract the first position
        disp = disp.^2; % square every cumulative displacement
    
        disp_xy = sum(disp,2); % sum x- and y- squared displacements
        disp_x = disp(:,1); % NOTE: squared x-displacements only
        disp_y = disp(:,2); % NOTE: squared y-displacements only
    
        dispMatXY(1:len-1,idx) = disp_xy;
        dispMatX(1:len-1,idx) = disp_x;
        dispMatY(1:len-1,idx) = disp_y;
    end
    
    MSD_xy = zeros(maxLen-1,1);
    MSD_x = zeros(maxLen-1,1);
    MSD_y = zeros(maxLen-1,1);
    for idx = 1:maxLen-1
        MSD_xy(idx) = sum(dispMatXY(idx,:))/nnz(dispMatXY(idx,:));
        MSD_x(idx) = sum(dispMatX(idx,:))/nnz(dispMatX(idx,:));
        MSD_y(idx) = sum(dispMatY(idx,:))/nnz(dispMatY(idx,:));
    end

end