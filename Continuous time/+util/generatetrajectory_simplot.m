%% Continuous-time trajectories with motion blur and localisation noise

% this function outputs positions rather than displacements, which is
% useful for e.g. plotting purposes

% INPUT:
% nSteps - number of displacements in each track
% tau - sampling time interval [ms]
% M - number of "subtimes" for each sampling time interval
% mC - threshold for the shutter function (should be <=M). =M for
% continuous illumination and open shutter mode.
% D1 - diffusion constant for state 1
% D2 - diffusion constant for state 2
% k12 - rate constant 1-->2 [1/ms]
% k21 - rate constant 2-->1 [1/ms]
% sigmaE - standard dev for the error in localization

% OUTPUT:
% out - a realistic 2-state particle trajectory in two dimensions

% Dependencies:
% shutter_func.m
% twoState_Markov.m

function [xRe,yRe,xMBLocNoise,yMBLocNoise,stateMat] =...
         generatetrajectory_simplot(nSteps, tau, M, D1, D2, k12, k21,sigmaE)
    
    import('util.shutter_func')
    import('util.twoState_Markov')

    % Derived variables
    nPos = nSteps+1; % number of positions
    % NOTE: unit conversion for tau below, from ms to s
    subSigmaB1 = sqrt(2*D1*tau)/sqrt(M); % std dev in displacements for each subtime interval
                                    % in particle state 1
    subSigmaB2 = sqrt(2*D2*tau)/sqrt(M); % std dev in displacements for each subtime interval
                                    % in particle state 2

    % --------------- Generate state sequences ------------------

    % Sampling-time sequence
%     p12 = k12/(k12+k21)*(1-exp(-(k12+k21)*tau)); % transition probability 1-->2
%     p21 = k21/(k12+k21)*(1-exp(-(k12+k21)*tau)); % transition probability 2-->1
%     stateVec = twoState_Markov(p12,p21,nPos);
    
    % Sub-time sequences
    dtau = tau/M; % sub-sampling time
    p12sub = k12/(k12+k21)*(1-exp(-(k12+k21)*dtau)); % sub-time transition probability 1-->2
    p21sub = k21/(k12+k21)*(1-exp(-(k12+k21)*dtau)); % sub-time transition probability 2-->1
    stateVec = util.twoState_Markov(p12sub,p21sub,nPos*M); % full, "unrolled" state sequence
    stateMat = reshape(stateVec',[M,nPos])';

    % -------- Generate Brownian displacements and positions ------------

    % Create binary masks for the two states
    mask1 = stateMat==1;
%     mask1 = stateVec==1; % logical indexing of state 1
%     mask1 = repmat(mask1', 1, M); % rows correspond to sampling times
%                                   % columns correspond to subtime intervals
    mask2 = ~mask1; % the corresponding indexing of state 2

    % Draw positions for the x-direction
    Dx1 = subSigmaB1*randn(nPos,M).*mask1; % draw displacements for state 1
    Dx2 = subSigmaB2*randn(nPos,M).*mask2; % draw displacements for state 2
    Dx = Dx1+Dx2;       % holds all actual displacements
                        % rows = different sampling times (labeled by n)
                        % columns = different subtime intervals (labeled by m)                
    DxRe = reshape(transpose(Dx),[1,nPos*M]); 
                        % reshape so that we get one long array of displacements
                        % (labeled by q = M*(n-1) + m)
    xRe = cumsum(DxRe); % actual positions for different q 
    x = transpose(reshape(xRe,[M,nPos])); % matrix with elements x_{n,m}

    % Draw positions for the y-direction
    Dy1 = subSigmaB1*randn(nPos,M).*mask1; % draw displacements for state 1
    Dy2 = subSigmaB2*randn(nPos,M).*mask2; % draw displacements for state 2
    Dy = Dy1+Dy2;       % holds all actual displacements
                        % rows = different sampling times (labeled by n)
                        % columns = different subtime intervals (labeled by m)                
    DyRe = reshape(transpose(Dy),[1,nPos*M]); 
                        % reshape so that we get one long array of displacements
                        % (labeled by q = M*(n-1) + m)
    yRe = cumsum(DyRe); % actual positions for different q 
    y = transpose(reshape(yRe,[M,nPos])); % matrix with elements x_{n,m}

    % Apply motion blur to Brownian positions 
    xMB = zeros(nPos,1);  % motion blurred x-positions
    yMB = zeros(nPos,1);  % motion blurred y-positions
    s = util.shutter_func(1:M,M);
    for n=1:nPos
        xMB(n) = sum(s.*x(n,:));
        yMB(n) = sum(s.*y(n,:));
    end

    % Apply localization noise
    xMBLocNoise = xMB + sigmaE*randn(nPos,1); % add localization noise to the x-positions
    yMBLocNoise = yMB + sigmaE*randn(nPos,1); % add localization noise to the y-positions

    % positions = [xMBLocNoise,yMBLocNoise]; % position data
    % Displacement which include both motion blur and localization noise 
    % Dx_final = x(2:nPos) - x(1:nPos-1); % x-displacement
    % Dy_final = y(2:nPos) - y(1:nPos-1); % y-displacement
    % 
    % displacements = sqrt(Dx_final.^2 + Dy_final.^2); % 2D displacement

    % new code:
    % x = x(:,1); % pick out positions at the beginning of each time interval
    % y = y(:,1); % pick out positions at the beginning of each time interval

end

