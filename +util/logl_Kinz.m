function enslogL = logl_Kinz(paramArr, exPars, data)
    % 
    % Approximative likelihood calculation for fast switching
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
    if D1 <= 0 || D2 <= 0 || D1 <= D2 || p12 <= 0 || p12 >= 1 || p21 <= 0 || p21 >= 1 % these inputs are not allowed (because they are unphysical)
        
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

            delta = data{i}; % current trajectory
 
            % define likelihood contribution functions
            % Lstay_func = @(y) pi_1*exp(-tau*k12)*1/(4*pi*tau*D1)*exp(-delta.^2/(4*tau*D1)) +... % begin and stay in state 1
            %                   pi_2*exp(-tau*k21)*1/(4*pi*tau*D2)*exp(-delta.^2/(4*tau*D2)); % begin and stay in state 2
            Lstay_func = @(y) pi_1*exp(-tau*k12)*1/(4*pi*tau*D1*(1-2*Rmb)+2*sigmaE^2)*exp(-delta.^2/(4*tau*D1*(1-2*Rmb)+2*sigmaE^2)) +... % begin and stay in state 1
                              pi_2*exp(-tau*k21)*1/(4*pi*tau*D2*(1-2*Rmb)+2*sigmaE^2)*exp(-delta.^2/(4*tau*D2*(1-2*Rmb)+2*sigmaE^2)); % begin and stay in state 2
            
            % Lswitch_func = @(w,y) tau*1/(4*pi*tau*(w*D1+(1-w)*D2)).*exp(-delta.^2/(4*tau*(w*D1+(1-w)*D2)))*exp(-tau*(k12*w+k21*(1-w)))*...
            %                       ((pi_1*k12*(1+tau*k21*w)+pi_2*k21*(1+tau*k12*(1-w)))*besseli(0,2*tau*sqrt(k12*k21*w.*(1-w)))-...
            %                       (tau*k12*k21)*(pi_1*w+pi_2*(1-w))*besseli(2,2*tau*sqrt(k12*k21.*w.*(1-w))));
            Lswitch_func = @(w,y) tau*1/(4*pi*tau*(w*D1+(1-w)*D2)*(1-2*Rmb)+2*sigmaE^2).*exp(-delta.^2/(4*tau*(w*D1+(1-w)*D2)*...
                                  (1-2*Rmb)+2*sigmaE^2))*exp(-tau*(k12*w+k21*(1-w)))*...
                                  ((pi_1*k12*(1+tau*k21*w)+pi_2*k21*(1+tau*k12*(1-w)))*besseli(0,2*tau*sqrt(k12*k21*w.*(1-w)))-...
                                  (tau*k12*k21)*(pi_1*w+pi_2*(1-w))*besseli(2,2*tau*sqrt(k12*k21.*w.*(1-w))));
            
            % evaluate likelihood
            Lstay = Lstay_func(delta);

            epsilon_abs = 1e-10; % absolute tolerance for numerical integration
            epsilon_rel = 1e-6; % relative tolerance for numerical integration
            Lswitch = integral(@(w) Lswitch_func(w,delta), 0,1, 'ArrayValued',true,...
                'AbsTol',epsilon_abs,'RelTol',epsilon_rel);
            
            logL = log(Lstay+Lswitch); % likelihood of each step
            logL = sum(logL); % likelihood of trajectory
    
        if isnan(logL) % error handling
            error('reaction rates are too high')
        end

        enslogL = enslogL + logL; % increment ensemble likelihood
        end        
    end
end


