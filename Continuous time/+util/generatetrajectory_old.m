%% Generate trajectories for a continuous-time 2-state process

% INPUT:
% k12 - rate constant 1 --> 2
% k21 - rate constant 2 --> 1
% D1 - diffusion constant for state 1
% D2 - diffusion constant for state 2
% tau - sampling time
% nSteps - number of sampling intervals = number of displacements

% OUTPUT:
% trackObs - a trajectory of displacements from 2-dimensional diffusion

% Dependencies:
% none

function trackObs = generatetrajectory(k12,k21,D1,D2,tau,nSteps)

    nIntervals = nSteps; % draw one random number for each time interval
    stopTime = nIntervals*tau;
    pi_1 = k21/(k12+k21); % stationary probability to be in state 1
%     pi_2 = 1-pi_1; % stationary probability to be in state 2

    % Generate times for when state transitions occur
    stateVec = [];
    timeVec = [0]; % times corresponding to stateVec, starting from time 0
    time = 0;
    u = rand;
    if u <= pi_1 % start in state 1
        
        stateVec(1) = 1;
        
    else % start in state 2
    
        stateVec(1) = 2;
    end
    
    while time <= stopTime
        
        if stateVec(end) == 1
            
%             time =  time + exprnd(1/k12);
            time =  time -1/k12*log(rand);
            timeVec(end+1) = time;
            stateVec(end+1) = 2;
            
        else 
        
%             time = time + exprnd(1/k21);
            time =  time -1/k21*log(rand);
            timeVec(end+1) = time;
            stateVec(end+1) = 1;
        end
    end
    timeVec(end) = stopTime; % truncate the last time to the stopping time
    stateVec(end) = []; % remove the last state occuring after stopTime

    % handle case with no state transitions
    if isempty(timeVec(2:end-1))

        if stateVec(1) == 1
        
            wVec = ones(1,nIntervals); % in state 1 during the entire trajectory
        
        else

            wVec = zeros(1,nIntervals); % in state 2 during the entire trajectory
        end

        % Draw displacements based on weight vector
        Dvec = wVec*D1 + (1-wVec)*D2; % weighted diffusion constants for each step
%       trackObs = sqrt(2*Dvec*tau).*raylrnd(1,1,nSteps); % corresponding displacements
        trackObs = raylrnd(sqrt(2*Dvec*tau));

    else % continue as usual

        % Handle sub-time transitions (this should be done before the next partitioning which assumes no sub-time transitions)
        tranInd = ceil(timeVec(2:end-1)/tau); % find transition intervals
        uniTranInd = unique(tranInd); % intervals of transitions without doublets,triplets etc.
        remInd = [];
        altInd = {};
        wVecN = zeros(1,nIntervals); % hold all weights stemming from subtransitions of an even number
        
        for k = 1:length(uniTranInd)
        
            tran = uniTranInd(k); % current unique transition interval
            y = sum(tranInd == tran); % number of sub-time transitions in this interval
            if y > 1 % sub-time transitions
        
                % compute all sub-times in this interval
                ti = find(tranInd==tran); % all indices corresponding to transitions within the same interval
                rt = zeros(1,y+1); % time differences between each sub-time transition and sampling interval
                rt(1) = timeVec(1+ti(1))-(tranInd(ti(1))-1)*tau;
                rt(2:end-1) = diff(timeVec(1+ti));
                rt(end) = (tranInd(ti(1))-1)*tau+tau - timeVec(1+ti(end));
        
                % compute the sub-time weights for the interval
                si = [stateVec(ti),stateVec(ti(end)+1)]; % corresponding states
                si(si==2) = 0;
                wVecI = si.*rt/tau;
                w1 = sum(wVecI); % weight of time in state 1 in this sumtime interval
                
                % store this information to later alter timeVec, stateVec and tranInd
                if stateVec(ti(1)) == 1 % previous state, before first sub-time transition
        
                    if tranInd(ti(1))-1 == 0 % sub-time transitions in the first interval
        
                        tNew = w1*tau; % position of new point
        
                    else % transition in any but the first interval
                    
                        tNew = (tranInd(ti(1))-1)*tau + w1*tau; % position of new point
                    end
        
                else % previous state = 2, before first sub-time transition
        
                    if tranInd(ti(1))-1 == 0 % sub-time transitions in the first interval
        
                        tNew = (1-w1)*tau;
        
                    else
                    
                        tNew = (tranInd(ti(1))-w1)*tau; % position of new point
                    end
                end
        
                if mod(y,2) == 1 % an odd number of transitions in this interval
                
                    remInd = [remInd,ti(2:end)]; % indices of values to be removed. Leave one index for later modification.
                    altInd{end+1} = [ti(1),tNew,1]; %  [index to modify, new value, odd/even]
                
                else % an even number of transitions in this interval
        
                    wVecN(tran) = w1; % add weight to this interval
                    remInd = [remInd,ti]; % indices of values to be removed. Remove all indices.
                    altInd{end+1} = [ti(1),tNew,0]; %  [index to modify, new value, odd/even]
                end
            end
        end
        
        % remove sub-time transitions and replace by a single transition yielding the same weights
        tranIndN = 1:length(tranInd);
%         bina = ones(1,length(tranInd)); % mask for removing extra sub-time transitions
%         for k = 1:length(remInd)
%         
%             if ismember(remInd(k),tranIndN)
%         
%                 bina(remInd(k)) = 0;
%             end
%         end
        bina = ~ismember(tranIndN,remInd);
        timeVecN = timeVec(2:end-1).*bina; % NOTE: shortening this vector now, so not to repeat that when partitioning into sampling times as usual
        
        for k = 1:length(altInd)
            
            arr = altInd{k};
        
            if arr(3) == 1 % an odd number of subtransitions
                
                timeVecN(arr(1)) = arr(2); % update value
            end % for intervals with an even number of transitions, we remove all transitions
        end
        timeVecN = timeVecN(timeVecN ~= 0); % remove extra transitions
        stateVecN = (stateVec(1:length(timeVecN)+1));

        % handle the case were an even number of transitions happen in a single time interval, and no other transitions occur
        if isempty(timeVecN)
            
            if stateVec(1) == 1
        
                wVec = ones(1,nIntervals); % in state 1 during the entire trajectory
        
            else

                wVec = zeros(1,nIntervals); % in state 2 during the entire trajectory
            end
        else % continue as usual
         
            % Partition into sampling times
            timeVec = timeVecN; % use the sub-transitioning result
            stateVec = stateVecN; % use the sub-transitioning result
            
            % timeVec = timeVec(2:end-1);
            tranInd = ceil(timeVec/tau); % find transition intervals
            restTimes1 = timeVec-((tranInd-1)*tau); % times after last sampling time
            restTimes2 = tau-restTimes1; % times before next sampling time
            
            wVecPrel = restTimes1/tau.*(stateVec(1:end-1) == 1) + restTimes2/tau.*(stateVec(1:end-1)==2); % vector of weights for state 1
            
            wVec = zeros(1,nIntervals); % final vector of weights
            for k = 1:length(wVecPrel) % add weights at the correct interval
            
                wVec(tranInd(k)) = wVecPrel(k);
            end
            
            % add weights in-between switches
            yVec = zeros(1,nIntervals); % binary state vector
            yVec(1:tranInd(1)) = 1*(stateVec(1)==1) + 0*(stateVec(1)==2);
            for k = 1:length(tranInd)-1 % add weights at the correct interval
            
                yVec(tranInd(k):tranInd(k+1)) = 1*(stateVec(k+1)==1) + 0*(stateVec(k+1)==2); % NOTE: 0-assign for state 2 is redundant
            end
            yVec(tranInd(length(tranInd)):end) = 1*(stateVec(end)==1) + 0*(stateVec(end)==2);
            
            wVec(wVec==0) = yVec(wVec==0);
            upN = find(wVecN ~= 0); % intervals with an even number of sub-time transitions
            for i = 1:length(upN)
            
                if wVec(i) == 1
            
                    wVec(i) = wVec(i)*wVecN(i);
                else
                    wVec(i) = wVec(i)+wVecN(i);
                end
            end
        end

        % Draw displacements based on weight vector
        Dvec = wVec*D1 + (1-wVec)*D2; % weighted diffusion constants for each step
%       trackObs = sqrt(2*Dvec*tau).*raylrnd(1,1,nSteps); % corresponding displacements
        trackObs = raylrnd(sqrt(2*Dvec*tau));
    end
end






