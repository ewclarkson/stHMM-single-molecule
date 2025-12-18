function enslogL = logl_Kinz_signal(paramArr, exPars, data)
    % 
    % Approximative likelihood calculation for fast switching for a
    % model in which the means of gaussian emission distributions are 
    % unkown as well as the state switching rates. 
    % 
    % assume that mu1 > mu2 and that sigma1 = sigma2.
    % 
    % Input:
    % mu1     - mean associated with state 1
    % mu2     - mean associated with state 2
    % k12     - transition rate from state 1 to state 2
    % k21     - transition rate from state 2 to state 1
    % exPars - cell array of experimental parameters.
    %          exPars = {'tau', 5; 'sigma1', 1; 'sigma2', 1}
    %          where 'tau' is the sampling time, 'Rmb' is the motion blur
    %          coefficient and 'sigmaE' is the localisation error std.
    %          NOTE: simply set Rmb=sigmaE=0 in case no nosie corrections
    %          are to be used.
    % data - cell array of signal measurements
    % 
    % Output:
    % loglikelihood - the exact logarithmic likelihood of the signal, 
    %                 given the input parameters 
    % 
    % Dependencies:
    % logsumexp2.m

    % extract experimental parameters
    for i = 1:size(exPars,1)
        specPar = exPars{i,1}; % name of specific parameter

        if strcmp(specPar,'tau')
            tau = exPars{i,2}; % sampling time

        elseif strcmp(specPar,'sigma')
            sigma = exPars{i,2}; % motion blur coefficient

        % elseif strcmp(specPar,'sigma2')
        %     sigma2 = exPars{i,2}; % localisation error std
        end
    end
    
    % sigma2 = sigma1; % otherwise, include broadening in likelihood function

    % extract parameters
    mu1 = paramArr(1);
    mu2 = paramArr(2);
    k12 = paramArr(3);
    k21 = paramArr(4);

    p12 = k12/(k12+k21)*(1-exp(-(k12+k21)*tau)); % transition probability 1 --> 2
    p21 = k21/(k12+k21)*(1-exp(-(k12+k21)*tau)); % transition probability 2 --> 1
%     if  D1 <= 0 || D2 <= 0 || D1 <= D2 || k12 <= 0 || k21 <= 0 % forbidden inputs
    if mu1 <= mu2 || p12 < 0 || p12 >= 1 || p21 < 0 || p21 >= 1 % these inputs are not allowed (because they are unphysical)
        
        enslogL = -Inf; % return a log-likelihood of -Inf
    
    % elseif k12*tau > 2 || k21*tau > 2 % a very high mean number of transitions
    % 
    %     enslogL = -Inf; % ignore very high rates
        
    else % continue with the algorithm

        enslogL = 0; % accumulated log-likelihood

        % --------------------- define constants and functions ------------------------

        pi_1 = k21/(k12+k21);
        pi_2 = k12/(k12+k21); % =1-pi1
    
        for i = 1:length(data) % do for every trajectory

            y = data{i}; % current trajectory
 
            % define likelihood contribution functions
            Lstay_func = @(y) pi_1*exp(-tau*k12)*1/sqrt(2*pi*sigma^2)*exp(-(y-mu1).^2/(2*sigma^2)) +... % begin and stay in state 1
                            pi_2*exp(-tau*k21)*1/sqrt(2*pi*sigma^2)*exp(-(y-mu2).^2/(2*sigma^2)); % begin and stay in state 2
            
            Lswitch_func = @(w,y) tau*1/sqrt(2*pi*sigma^2).*exp(-(y-(w*mu1+(1-w)*mu2)).^2./(2*sigma^2))*exp(-tau*(k12*w+k21*(1-w)))*...
                                  ((pi_1*k12*(1+tau*k21*w)+pi_2*k21*(1+tau*k12*(1-w)))*besseli(0,2*tau*sqrt(k12*k21*w.*(1-w)))-...
                                  (tau*k12*k21)*(pi_1*w+pi_2*(1-w))*besseli(2,2*tau*sqrt(k12*k21.*w.*(1-w))));
       
            % evaluate likelihood
            Lstay = Lstay_func(y);

            epsilon_abs = 1e-10; % absolute tolerance for numerical integration
            epsilon_rel = 1e-6; % relative tolerance for numerical integration
            Lswitch = integral(@(w) Lswitch_func(w,y), 0,1, 'ArrayValued',true,...
                'AbsTol',epsilon_abs,'RelTol',epsilon_rel);
            
            logL = log(Lstay+Lswitch); % likelihood of each step
            logL = sum(logL); % likelihood of trajectory

        if isnan(logL) % error handling
            error('reaction rates are too high!')
        end

        enslogL = enslogL + logL; % increment ensemble likelihood
        end        
    end
end


