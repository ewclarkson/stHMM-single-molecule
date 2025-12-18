function [x,y] = traj1state_pos(nPos, tau, D, sigmaE)
    % 
    % Two-dimensional diffusion with localisation errors
    % 
    % Input:
    % nPos = number of positions
    % tau = sampling time
    % D = diffusion constant
    % sigmaE = localisation error std
    % 
    % Output:
    % x = x-coordinates of random walk
    % y = y-coordinates of random walk
    % 

    sigma = sqrt(2*tau*D);
    x = cumsum(sigma*randn(nPos,1) + sigmaE*randn(nPos,1));
    y = cumsum(sigma*randn(nPos,1) + sigmaE*randn(nPos,1));

end