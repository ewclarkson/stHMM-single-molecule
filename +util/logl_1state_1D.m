function enslogL = logl_1state_1D(paramArr, exPars, data)
    % 
    % Discrete-time likelihood calculation with a forward algorithm.
    % For 1-dimensional diffusion!
    % 
    % Input:
    % paramArr - [D1] diffusion constant associated with state 1
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

    if D1 <= 0 % forbidden input, which is unphysical (should be taken care of by the priors)
        
        enslogL = -Inf; % return a log-likelihood of -Inf
        
    else % continue with the likelihood

        enslogL = 0;
        logL_func = @(disp, D, deltaT) -0.5*log(4*pi*D*deltaT*(1-2*Rmb)+2*sigmaE^2) - ...
                                            disp.^2/(4*D*deltaT*(1-2*Rmb)+2*sigmaE^2);
        for i = 1:length(data) % do for every trajectory
    
            trackDisps = data{i}; % current trajectory
            enslogL = enslogL + sum(logL_func(trackDisps, D1, tau));
        end
    end

