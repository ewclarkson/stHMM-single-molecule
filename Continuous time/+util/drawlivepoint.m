function [uposNew, logLnew, posNew, reject] = drawlivepoint(livePoints, logLworst, exPars, data, priorPars, loglfun)

    % 
    % Uniformly draw a point obeying a likelihood constraint, by estimating 
    % an ellipsoidal isocontour bounding the live points. Accept a proposed 
    % point if it has a likelihood higher than the constraint. This method 
    % is presented in Shaw et al "Efficient Bayesian inference for 
    % multimodal problems in cosmology" Mon. Not. R. Astron. Soc. 378, 
    % 1365–1370 (2007).
    % 
    % Input:
    % livePoints - structure array of all live points with one field 'pos'
    %              that contains all positions in parameter space.
    %              livePoints(i).pos = [D1_i, D2_i, k12_i, k21_i] for point i                
    % logLworst  - log-likelihood constraint, a lower bound
    % exPars     - cell array of experimental parameters.
    %              exPars = {'tau', 5; 'Rmb', 1/6; 'sigmaE', 0.1}
    %              where 'tau' is the sampling time, 'Rmb' is the motion blur
    %              coefficient and 'sigmaE' is the localisation error std.
    %              NOTE: simply set Rmb=sigmaE=0 in case no nosie corrections
    %              are to be used.
    % data       - cell of trajectories of displacements
    % prior      - prior distribution parameters
    % loglfun    - handle to likelihood function. @chosenlikelihoodfunction
    % 
    % Output:
    % posNew  - position of new point
    % logLnew - likelihood of new point
    % reject  - number of rejections of proposed points
    % 
    % Dependencies:
    % logsumexp2.m

    % Preliminaries
    pointArr = reshape([livePoints.upos],4,[])'; % numerical array of all live points' positions
    [nPoints, nDims] = size(pointArr); % extract #points and #dimensions
    
    % ------------- approximate a bounding ellipsoid ---------------------
    
    mean_p = mean(pointArr); % estimated mean of points
    
    % Estimate covariance matrix statistically
    app_C = cov(pointArr,1); % compute covariance matrix
    inv_appC = inv(app_C); % inverse of covariance matrix
    
    % Find a bounding ellipsoid
    kVec = zeros(1,nPoints); % helps define an ellipsoid
    for i = 1:nPoints % compute a value in kVec for every point
    
        kVec(i) = (pointArr(i,:)-mean_p)*inv_appC*(pointArr(i,:)-mean_p)';
%         kVec(i) = (pointArr(i,:)-mean_p)*(app_C\((pointArr(i,:)-mean_p)'));
    end
    kScale = max(kVec); % pick the largest k
    kNew = kScale*1.06^2; % slightly enlarge the ellipsoid
    
    % ------------ sample uniformly from bounding ellipsoid --------------
    
    % Compute a transformation matrices from ball to ellipsoid
    [R,D] = eig(app_C); % matrix of eigenvectors and of eigenvalues
%     [~,D,R] = eig(app_C); 
%     R = R';
%     T = sqrt(kNew)*R'*sqrt(D)*R; % transformation matrix from ball to ellipsoid
    T = sqrt(kNew)*R*sqrt(D)*R'; % transformation matrix from ball to ellipsoid

    % -- test correctness --
    % zArr = zeros(nPoints,nDims);
    % for k = 1:nPoints % sample within D-ball
    %     zt = randn(1,nDims); % random vector of standard gaussians
    %     zt = (rand)^(1/nDims)*zt/norm(zt); % sample within unit ball
    %     zArr(k,:) = zt; % store z
    % end
    % 
    % yVec = zeros(nPoints,nDims);
    % for k = 1:nPoints
    %     yVec(k,:) = (T*zArr(k,:)')'+mean_p; % final random deviates
    % end
    % -- end of test --

    logLnew = logLworst-1;
    reject = -1; % number of rejections
    while logLnew <= logLworst % do until we are inside the region of L>Lworst

        % Generate a random deviate within D-ball
        z = randn(1,nDims); % random vector of standard gaussians
        z = rand^(1/nDims)*z/norm(z); % sample within unit ball
                
        % Apply transformation to D-ball random deviates
        u = (T*z')'+mean_p; % final random deviates
        if any(u<=0) || any(u>=1) % forbidden values
            logLnew = logLworst-1;
            reject = reject+1;
        else % continue as usual
            % Evaluate likelihood in original parameter space
            y = zeros(1,4);
    
            for idx = 1:4
                % y(idx) = unifinv(u(idx),priorPars(idx,1),priorPars(idx,2));
                if strcmp(priorPars{idx,1},'uniform')
                    y(idx) = unifinv(u(idx),priorPars{idx,2},priorPars{idx,3});
                elseif strcmp(priorPars{idx,1},'lognormal')
                    y(idx) = logninv(u(idx),priorPars{idx,2},priorPars{idx,3});
                end                
            end

            if y(2) > y(1)
                    y = [y(2),y(1),y(4),y(3)]; % enforce D1>D2
                    % ------- new code -------------
                    for k = 1:4
                        if strcmp(priorPars{k,1},'uniform')
                            u(k) = unifcdf(y(k),priorPars{k,2},priorPars{k,3});
                        
                        elseif strcmp(priorPars(k,1),'lognormal')
                            u(k) = logncdf(y(k),priorPars{k,2},priorPars{k,3});
                        end
                    end 
                    % --------- end of new code ------------- 
                    % u = [u(2),u(1),u(4),u(3)]; % enforce y(1)>y(2)
            end
    
            logLnew = loglfun(y, exPars, data); % likelihood of new point
            reject = reject+1;
        end
    end
    
    uposNew = u;
    posNew = y; % output
end

% scatter(pointArr(:,1), pointArr(:,2)) 
% hold on
% scatter(yVec(:,1), yVec(:,2))

% scatter(pointArr(:,3), pointArr(:,4))
% hold on
% scatter(yVec(:,3), yVec(:,4))

% scatter(pointArr(:,1), pointArr(:,3))
% hold on
% scatter(yVec(:,1), yVec(:,3))

% scatter(pointArr(:,1), pointArr(:,4))
% hold on
% scatter(yVec(:,1), yVec(:,4))

% scatter(pointArr(:,2), pointArr(:,3))
% hold on
% scatter(yVec(:,2), yVec(:,3))

% scatter(pointArr(:,2), pointArr(:,4))
% hold on
% scatter(yVec(:,2), yVec(:,4))

