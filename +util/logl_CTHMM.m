function enslogL = logl_CTHMM(paramArr, exPars, data)

    % use a modified version that handles log-likelihoods. Assume 2D tracks.

    % Input:
    % D1 - diffusion constant associated with state 1
    % D1 - diffusion constant associated with state 1
    % k12 - transition rate from state 1 to state 2
    % k21 - transition rate from state 2 to state 1
    % exPars - cell array of experimental parameters.
    %          exPars = {'tau', 5; 'Rmb', 1/6; 'sigmaE', 0.1}
    %          where 'tau' is the sampling time, 'Rmb' is the motion blur
    %          coefficient and 'sigmaE' is the localisation error std.
    %          NOTE: simply set Rmb=sigmaE=0 in case no nosie corrections
    %          are to be used.
    % data - cell array of track displacements
    
    % Output:
    % loglikelihood - the exact logarithmic likelihood of the track, given the input parameters 
    
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
    if D1 <= 0 || D2 <= 0 || D1 <= D2 || p12 < 0 || p12 >= 1 || p21 < 0 || p21 >= 1 % these inputs are forbidden
       
        enslogL = -Inf; % return a log-likelihood of -Inf
    
    % elseif k12*tau > 2 || k21*tau > 2
    % 
    %     enslogL = -Inf; % forbid very high rates
            
    else % continue with the algorithm

        enslogL = 0; % accumulated log-likelihood

        % --------------------- define constants ------------------------------
        
        % p12 = 1-exp(-k12*tau); % transition probability 1 --> 2 
        % p21 = 1-exp(-k21*tau); % transition probability 2 --> 1
        % log_p11 = log(1-p12); % log transition probability 1 --> 1 with no switches
        % log_p22 = log(1-p21); % log transition probability 2 --> 2 with no switches

        log_pi1 = log(k21/(k12+k21)); % stationary probability to be in state 1
        log_pi2 = log(k12/(k12+k21)); % stationary probability to be in state 2

        for j = 1:length(data)

            trackDisps = data{j}; % current trajectory
            N = length(trackDisps); % number of displacements in the track

                % ----------- define functions for later computation -----------
        
            % for 1 --> 1, with at least two switches
            L_switch_same_1 = @(w,delta) 1./(4*pi*tau*(w*D1+(1-w)*D2)*(1-2*Rmb)+2*sigmaE^2).*exp(-delta.^2./(4*tau*(w*D1+(1-w)*D2)*(1-2*Rmb)+2*sigmaE^2))*...
                tau.*exp(-tau*(k12*w+k21*(1-w)))*tau*k12*k21.*w.*...
                (besseli(0,2*tau*sqrt(k12*k21*w.*(1-w))) - besseli(2,2*tau*sqrt(k12*k21.*w.*(1-w))));
    % 
            % for 2 --> 2, with at least two switches
            L_switch_same_2 = @(w,delta) 1./(4*pi*tau*(w*D1+(1-w)*D2)*(1-2*Rmb)+2*sigmaE^2).*exp(-delta.^2./(4*tau*(w*D1+(1-w)*D2)*(1-2*Rmb)+2*sigmaE^2))*...
                tau.*exp(-tau*(k12*w+k21*(1-w)))*tau*k12*k21.*(1-w).*...
                (besseli(0,2*tau*sqrt(k12*k21*w.*(1-w))) - besseli(2,2*tau*sqrt(k12*k21*w.*(1-w))));
    % 
            % for at least one switch. Common for 1--> 2 and 2 --> 1, just add factors of k12 or k21 in the end
            L_switch_opp = @(w,delta) 1./(4*pi*tau*(w*D1+(1-w)*D2)*(1-2*Rmb)+2*sigmaE^2).*exp(-delta.^2./(4*tau*(w*D1+(1-w)*D2)*(1-2*Rmb)+2*sigmaE^2))*...
            tau.*exp(-tau*(k12*w+k21*(1-w))).*besseli(0,2*tau*sqrt(k12*k21*w.*(1-w)));
    
            % L_switch_opp_2 = @(w,delta) 1./(4*pi*tau*(w*D1+(1-w)*D2)).*exp(-delta.^2./(4*tau*(w*D1+(1-w)*D2)))*...
            %     tau.*exp(-tau*(k12*w+k21*(1-w)))*k12.*besseli(0,2*tau*sqrt(k12*k21*w.*(1-w)));
    %   
            % for 2 --> 1, with at least one switch
            % L_switch_opp_1 = @(w,delta) 1./(4*pi*tau*(w*D1+(1-w)*D2)).*exp(-delta.^2./(4*tau*(w*D1+(1-w)*D2)))*...
            %     tau.*exp(-tau*(k12*w+k21*(1-w)))*k21.*besseli(0,2*tau*sqrt(k12*k21*w.*(1-w)));
    
            logLstay_1 = -tau*k12-log(4*pi*D1*tau*(1-2*Rmb)+2*sigmaE^2)-trackDisps.^2/(4*D1*tau*(1-2*Rmb)+2*sigmaE^2); % for 1 --> 1
    
            logLstay_2 = -tau*k21-log(4*pi*D2*tau*(1-2*Rmb)+2*sigmaE^2)-trackDisps.^2/(4*D2*tau*(1-2*Rmb)+2*sigmaE^2); % for 2 --> 2
            
            epsilon_abs = 1e-10; % absolute tolerance for numerical integration
            epsilon_rel = 1e-6; % relative tolerance for numerical integration
            % for 1 --> 1, with at least two switches
            logLsame_w1 = integral(@(w) L_switch_same_1(w,trackDisps), 0,1, 'ArrayValued',true,...
                'AbsTol',epsilon_abs,'RelTol',epsilon_rel);
            logLsame_w1 = log(logLsame_w1);
    
            % for 2 --> 2, with at least two switches
            logLsame_w2 = integral(@(w) L_switch_same_2(w,trackDisps), 0,1, 'ArrayValued',true,...
                'AbsTol',epsilon_abs,'RelTol',epsilon_rel);
            logLsame_w2 = log(logLsame_w2);
    
            % for at least one switch
            logLopp = integral(@(w) L_switch_opp(w,trackDisps), 0,1, 'ArrayValued',true,...
                'AbsTol',epsilon_abs,'RelTol',epsilon_rel);
            logLopp_w2 = log(logLopp*k12);
            logLopp_w1 = log(logLopp*k21);
            
    
            % -------------------- compute all forward variables --------------------------
    
            % initialise forward variables
            alpha_1s = util.logsumexp2(log_pi1 + util.logsumexp2(logLstay_1(1),logLsame_w1(1)), log_pi2 + logLopp_w1(1)); % alpha1(1)
            alpha_2s = util.logsumexp2(log_pi2 + util.logsumexp2(logLstay_2(1),logLsame_w2(1)), log_pi1 + logLopp_w2(1)); % alpha1(2)
            alphaVec_1 = zeros(1,N); alphaVec_1(1) = alpha_1s; % initilise array for all values of alpha_1
            alphaVec_2 = zeros(1,N); alphaVec_2(1) = alpha_2s; % initilise array for all values of alpha_2
    
            for n = 2:N % do for the remaining displacements
                
                alphaVec_1(n) = util.logsumexp2(alphaVec_1(n-1) + util.logsumexp2(logLstay_1(n),logLsame_w1(n)), alphaVec_2(n-1)+logLopp_w1(n)); 
                alphaVec_2(n) = util.logsumexp2(alphaVec_2(n-1) + util.logsumexp2(logLstay_2(n),logLsame_w2(n)), alphaVec_1(n-1)+logLopp_w2(n)); 
            end    
    
            logL = util.logsumexp2(alphaVec_1(N), alphaVec_2(N)); % total log-likelihood
    
            % error handling
            if isnan(logL)
        
                error('reaction rates are too high')
            end
    
            enslogL = enslogL + logL; % increment ensemble likelihood
        end
    end 
end



