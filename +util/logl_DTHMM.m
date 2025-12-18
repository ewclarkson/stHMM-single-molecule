function enslogL = logl_DTHMM(paramArr, exPars, data)
    % 
    % Discrete-time likelihood calculation with a forward algorithm
    % 
    % Input:
    % D1     - diffusion constant associated with state 1
    % D1     - diffusion constant associated with state 1
    % k12    - transition rate from state 1 to state 2
    % k21 - transition rate from state 2 to state 1
    % exPars - cell array of experimental parameters.
    %          exPars = {'tau', 5; 'Rmb', 1/6; 'sigmaE', 0.1}
    %          where 'tau' is the sampling time, 'Rmb' is the motion blur
    %          coefficient and 'sigmaE' is the localisation error std.
    %          NOTE: simply set Rmb=sigmaE=0 in case no nosie corrections
    %          are to be used.
    % data - cell array of track displacements
    % 
    % Output:
    % loglikelihood - the exact logarithmic likelihood of the track, given the input parameters 
    % 
    % Dependencies:
    % logsumexp2.m

    % extract experimental parameters
    for i = 1:size(exPars,1)
        specPar = exPars{i,1}; % name of specific parameter

        if strcmp(specPar,'tau')
            tau = exPars{i,2}; % sampling time

        elseif strcmp(specPar,'Rmb')
            Rmb = exPars{i,2}; % motion blur coefficient

        elseif strcmp(specPar,'sigmaE')
            sigmaE = exPars{i,2}; % localisation error std
        end
    end
    
    % extract parameters
    D1 = paramArr(1);
    D2 = paramArr(2);
    k12 = paramArr(3);
    k21 = paramArr(4);

    p12 = k12/(k12+k21)*(1-exp(-(k12+k21)*tau)); % transition probability 1 --> 2
    p21 = k21/(k12+k21)*(1-exp(-(k12+k21)*tau)); % transition probability 2 --> 1
%     if  D1 <= 0 || D2 <= 0 || D1 <= D2 || k12 <= 0 || k21 <= 0 % forbidden inputs
    if D1 <= 0 || D2 <= 0 || D1 <= D2 || p12 < 0 || p12 >= 1 || p21 < 0 || p21 >= 1 % these inputs are not allowed (because they are unphysical)
        
        enslogL = -Inf; % return a log-likelihood of -Inf
    
    % elseif k12*tau > 2 || k21*tau > 2 % a very high mean number of transitions
    % 
    %     enslogL = -Inf; % ignore very high rates
        
    else % continue with the algorithm

        enslogL = 0; % accumulated log-likelihood

        % --------------------- define constants and functions ------------------------

        log_p12 = log(p12); % log
        log_p21 = log(p21); % log
        log_p11 = log(1-p12); % log 
        log_p22 = log(1-p21); % log

        log_pi1 = log(p21/(p12+p21)); % the fraction of steps that the particle is in state 1
        log_pi2 = log(p12/(p12+p21)); % the fraction of steps that the particle is in state 2

        % logL_func = @(disp_j, D_i, deltaT) -log(4*pi*D_i*deltaT)-disp_j.^2/(4*D_i*deltaT); % the likelihood for every individual displacement
        logL_func = @(disp_j, D_i, deltaT) -log(4*pi*D_i*deltaT*(1-2*Rmb)+2*sigmaE^2) - ...
                                            disp_j.^2/(4*D_i*deltaT*(1-2*Rmb)+2*sigmaE^2); 
        % the likelihood for every individual displacement
    
        for i = 1:length(data) % do for every trajectory

            trackDisps = data{i}; % current trajectory
            N = length(trackDisps); % number of displacements

            % ----------- compute log-likelihoods -----------
    
            likeli_1 = logL_func(trackDisps, D1, tau); % emission probabilities all track displacements, given D1
            likeli_2 = logL_func(trackDisps, D2, tau); % emission probabilities all track displacements, given D2
      
            % -------------------- compute all alphas's --------------------------
    
            alpha_11 = log_pi1 + likeli_1(1); % compute alpha1(1)
            alpha_12 = log_pi2 + likeli_2(1); % compute alpha1(2)
            alphaArr_1 = zeros(1,N); alphaArr_1(1) = alpha_11; % initilise array for all values of alpha_1
            alphaArr_2 = zeros(1,N); alphaArr_2(1) = alpha_12; % initilise array for all values of alpha_2
    
            for j = 2:N % do for step 2 to step N
                
                alphaArr_1(j) = util.logsumexp2(alphaArr_1(j-1) + log_p11, alphaArr_2(j-1) + log_p21) + likeli_1(j); 
                % the prob. of observing displacement sequence 1,2,..j and % being in state 1 at step j given the model parameters

                alphaArr_2(j) = util.logsumexp2(alphaArr_1(j-1) + log_p12, alphaArr_2(j-1) + log_p22) + likeli_2(j); 
                % the prob. of observing displacement sequence 1,2,..j and being in state 2 at step j given the model parameters
            end    
    
            logL = util.logsumexp2(alphaArr_1(N), alphaArr_2(N));
    
        if isnan(logL) % error handling
            error('reaction rates are too high')
        end

        enslogL = enslogL + logL; % increment ensemble likelihood
        end        
    end
end


