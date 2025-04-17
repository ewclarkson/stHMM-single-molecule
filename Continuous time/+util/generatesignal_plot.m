%% Continuous-time signals from a 2-state single molecule experiment



function [sigOut,stateMat] = generatesignal_plot(nData, tau, M, mC, mu1, mu2, sigma1, sigma2, k12, k21)
    
    import('util.shutter_func')
    import('util.twoState_Markov')

    % Derived variables
    % nPos = nData+1; % number of positions
    % NOTE: unit conversion for tau below, from ms to s
    subSigmaB1 = sigma1/sqrt(M); % std dev in displacements for each subtime interval
                                    % in particle state 1
    subSigmaB2 = sigma2/sqrt(M); % std dev in displacements for each subtime interval
                                    % in particle state 2

    % --------------- Generate state sequences ------------------

    % Sub-time sequences
    dtau = tau/M; % sub-sampling time
    p12sub = k12/(k12+k21)*(1-exp(-(k12+k21)*dtau)); % sub-time transition probability 1-->2
    p21sub = k21/(k12+k21)*(1-exp(-(k12+k21)*dtau)); % sub-time transition probability 2-->1
    stateVec = util.twoState_Markov(p12sub,p21sub,nData*M); % full, "unrolled" state sequence
    stateMat = reshape(stateVec',[M,nData])';



    % -------- Generate Brownian displacements and positions ------------

    % Create binary masks for the two states
    mask1 = stateMat==1;
%     mask1 = stateVec==1; % logical indexing of state 1
%     mask1 = repmat(mask1', 1, M); % rows correspond to sampling times
%                                   % columns correspond to subtime intervals
    mask2 = ~mask1; % the corresponding indexing of state 2

    % Draw positions for the x-direction
    Dsig1 = (mu1/M+subSigmaB1*randn(nData,M)).*mask1; % draw signals for state 1
    Dsig2 = (mu2/M+subSigmaB2*randn(nData,M)).*mask2; % draw signals for state 2
    Dsig = Dsig1 + Dsig2; % all signal values
    sigOut = sum(Dsig,2)';
    
    % DsigRe = reshape(transpose(Dsig),[1,nData*M]); 
    %                     % reshape so that we get one long array of displacements
    %                     % (labeled by q = M*(n-1) + m)
    % sigRe = cumsum(DsigRe); % actual positions for different q 
    % sigOut = transpose(reshape(sigRe,[M,nData])); % matrix with elements x_{n,m}

    % Apply motion blur to positions 
    % xMB = zeros(nPos,1);  % motion blurred x-positions
    % yMB = zeros(nPos,1);  % motion blurred y-positions
    % s = shutter_func(1:M,mC);
    % for n=1:nPos
    %     xMB(n) = sum(s.*x(n,:));
    %     yMB(n) = sum(s.*y(n,:));
    % end
    % 

    % positions = x; % position data
     
end

