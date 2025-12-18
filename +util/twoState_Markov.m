%% Generate a discrete Markov chain

% that specifies the particle state (free or bound) at each time point. The
% initial state of the particle is chosen randomly according to the
% stationary probabilities of the two states
 
% INPUT:
% p12 - transition probability from state 1 to state 2
% p21 - transition probability from state 2 to state 1
% N - number of time points (no. of states to be generated)

% OUTPUT:
% out - a row vector containing all successive particle states at each time
% point. E.g. sVec = [1,1,2,2,2,1].


function out = twoState_Markov(p12, p21, N)

    pi_1 = p21/(p12+p21); % initialise stationary probabilities

    sVec = zeros(1,N); % initialise vector of states


    % --------------- assign the first state ---------------------------------

    u = rand(); % draw a rectangular random number

    if u <= pi_1
        sVec(1) = 1; % assign assign the first state (i.e. s1) to be state 1
    else
        sVec(1) = 2; % assign assign the first state (i.e. s1) to be state 2
    end


    % ----------------- assign all other states ------------------------------

    for i = 2:N
        u = rand(); % draw a rectangular random number

        if  sVec(i-1) == 1 && u <= p12 % conditional prob. from state 1 and transition prob. to state 2 are both fulfilled
            sVec(i) = 2; % make the transition from state 1 to state 2 at time step i

        elseif  sVec(i-1) == 2 && u <= p21 % conditional prob. from state 2 and transition prob. to state 1 are both fulfilled
            sVec(i) = 1; % make the transition from state 2 to state 1 at time step i

        else % the particle "chose" not to transition, i.e. to stay
            sVec(i) = sVec(i-1); % stay at the latest state at time step i
        end
    end  
  
  out = sVec;
  
end
    


