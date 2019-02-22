function state = matRad_doFSM(x,deltaT,options)

% from doi:10.1088/0031-9155/49/23/012

% state: number
% inhale    : 1
% exhale    : 2
% EOE       : 3
% irregular : 4

%% extract parameters
% number of points in decision line
K = round(options.timeLd./deltaT);
% state velocity ranges
velocityRangeS = [  options.velocityRangeInh;
                    options.velocityRangeExh;
                    options.velocityRangeEOE];
% threshold for velocity difference
cThetaVec = options.cTheta.*ones(4,1);
% threshold for amplitude
cLambdaVec = options.cLambda.*ones(8,1);
% threhold for state length
cSLength = options.cSLength;
% mapping function for the next state
fSNextMap = [2; 3; 1];
% time vector
t = (0:(numel(x)-1))'.*deltaT;

% vector containing most recent EOE velocities
vEOEVec = zeros(4,1);

%% starting values
% index of beginning of most recent lines
fmIndStart = ones(size(x));
% time of beginning of most recent lines
fmTStart = deltaT.*ones(size(x));
% knot position as a function of index (NOT MOST RECENT)
xKnot = zeros(size(x));
% breathing state of most recent lines
% assume we start in an exhale state
fmS = 2.*ones(size(x));

%% run state model
state = zeros(size(x));

% loop over entire stream
for i = K:numel(x)
    
    %% dynamic adaptation
    
    % calculate mean of the threshold vectors tracking the 4 most recent
    % periods
    cTheta = mean(cThetaVec);
    cLambda = mean(cLambdaVec);
    
    %% decision line fitting
    
    % get decision line index
    dIndStart = i-K+1;
    
    % do the fit
    dV = matRad_FSMgetFitVelocity(x(dIndStart:i),deltaT);
    
    if fmS(1) == 4
        % we are in an irregular breathing state
        
        % check if velocity of decision line is within range of either the
        % exhale or inhale state
        
        if ((velocityRangeS(1,1) < dV && dV < velocityRangeS(1,2)) || velocityRangeS(2,1) < dV && dV < velocityRangeS(2,2)) && i >= fmIndStart(1)+3
            
            %% enter regular breathing
            
            % now find out which state it is
            
            % first inhale state
            if velocityRangeS(1,1) < dV && dV < velocityRangeS(1,2)
                
                % set next state
                fSNext = 1;
            end
            
            % then exhale state
            if velocityRangeS(2,1) < dV && dV < velocityRangeS(2,2)
                
                % set next state
                fSNext = 2;
            end
            
            %% new line generation
            
            [fmTStart,fmIndStart,xKnot,fmS] = matRad_FSMnewLineGeneration(x,t,deltaT,fmTStart,fmIndStart,fmS,xKnot,i,K,fV,dV,fSNext);
            
            %% fitting optimization and dynamic adaptation
            
            [fmTStart,fmIndStart,xKnot,state,cLambdaVec,cThetaVec,vEOEVec] = matRad_FSMfittingOptimization(x,t,deltaT,fmTStart,fmIndStart,fmS,xKnot,state,cLambdaVec,cThetaVec,vEOEVec);
            
        end
        
        % note that we don't need to do anything extra if the velocity is
        % not in either of the inhale or exhale states (update without new
        % line segment)
        
    else
        % we are not in an irregular breathing state
        
        if i < fmIndStart(1)+2
            % we need at least two points to get a velocity in the most
            % recent line
            continue
        end
        
        % most recent line fitting
        fV = matRad_FSMgetFitVelocity(x(fmIndStart(1):(i-1)),deltaT);
        
        % calculate difference of decision and most recent line velocities
        fdTheta = abs(fV-dV);
        
        % check if velocity difference is greater than threshold
        if fdTheta > cTheta
            % the velocity difference is greater than the threshold; we may
            % be entering the next state
            % note that we don't need to do anything extra if the velocity
            % difference is less than the threshold (update without new
            % line segment)
            
            % check if velocity of decision line is within range of next
            % state
            fSNext = fSNextMap(fmS(1));
            if velocityRangeS(fSNext,1) < dV && dV < velocityRangeS(fSNext,2)
                % velocity is within range
                
                % calculate mean amplitude for the two most recent lines
                
                % checking if this is our first amplitude that we are
                % calculating
                if fmIndStart(1) == 1 || fmIndStart(2) >= fmIndStart(1)-1
                    % if it is, just reuse the most recent line
                    fdLambda = matRad_FSMgetAmplitude(x(fmIndStart(1):(i-1)),x(fmIndStart(1):(i-1)),deltaT);
                else
                    fdLambda = matRad_FSMgetAmplitude(x(fmIndStart(1):(i-1)),x(fmIndStart(2):fmIndStart(1)-1),deltaT);
                end
                
                % check if mean amplitude is greater than threshold
                if fdLambda > cLambda && i >= fmIndStart(1)+3
                    % greater than threshold
                    
                    %% new line generation
                    
                    [fmTStart,fmIndStart,xKnot,fmS] = matRad_FSMnewLineGeneration(x,t,deltaT,fmTStart,fmIndStart,fmS,xKnot,i,K,fV,dV,fSNext);
                    
                    %% fitting optimization and dynamic adaptation
                    
                    [fmTStart,fmIndStart,xKnot,state,cLambdaVec,cThetaVec,vEOEVec] = matRad_FSMfittingOptimization(x,t,deltaT,fmTStart,fmIndStart,fmS,xKnot,state,cLambdaVec,cThetaVec,vEOEVec);
                    
                else
                    % smaller than threhold
                    
                    %% local pruning
                    
                    [fmTStart,fmIndStart,fmS,xKnot] = matRad_FSMlocalPruning(fmTStart,fmIndStart,fmS,xKnot,velocityRangeS,dV);
                    
                end
                
            else
                % velocity is not within range
                
                % check if velocity of most recent line is within range of current state
                if ~(velocityRangeS(fmS(1),1) < fV && fV < velocityRangeS(fmS(1),2))
                    % velocity is not within range
                    % note that we don't need to do antyhing extra if the
                    % velocity is within range (update without new line
                    % segment)
                    
                    % check if velocity of decision line is within range of
                    % the inhale state
                    if velocityRangeS(1,1) < dV && dV < velocityRangeS(1,2)
                        % velocity is within range
                        
                        %% EOE enforcement
                        
                        [fmTStart,fmIndStart,xKnot,fmS,state,cLambdaVec,cThetaVec,vEOEVec] = matRad_FSMeoeEnforcement(x,t,deltaT,fmTStart,fmIndStart,fmS,xKnot,state,cLambdaVec,cThetaVec,vEOEVec,i,K);
                        
                    else
                        % velocity is not within range
                        
                        % check if velocity of decision line is within
                        % range of the exhale state
                        if velocityRangeS(2,1) < dV && dV < velocityRangeS(2,2)
                            % velocity is within range
                            
                            %% local pruning
                            
                            [fmTStart,fmIndStart,fmS,xKnot] = matRad_FSMlocalPruning(fmTStart,fmIndStart,fmS,xKnot,velocityRangeS,dV);
                            
                        else
                            % velocity is not within range
                            
                            %% enter irregular breathing
                            
                            % change state of most recent line
                            fmS(1) = 4;
                            
                            % the previous state might have also bee IRR, check to see and
                            % fix if so
                            if fmS(1) == fmS(2) && fmIndStart(2) > 1
                                
                                % update beginning times, indices, statesof the most
                                % recent lines, starting with the most recent one
                                % also knots
                                fmTStart    = circshift(fmTStart,-1);
                                fmIndStart  = circshift(fmIndStart,-1);
                                fmS         = circshift(fmS,-1);
                                xKnot       = circshift(xKnot,-1);
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
        % only do this check if not already in IRR
        if fmS(1) ~= 4 && i < fmIndStart(1)+2
            
            % redo most recent line fitting, but now including most recent data
            % point
            fV = matRad_FSMgetFitVelocity(x(fmIndStart(1):i),deltaT);
            
            % check if velocity of most recent line is within range of current state
            if ~(velocityRangeS(fmS(1),1) < fV && fV < velocityRangeS(fmS(1),2))
                % velocity is not within range
                
                %% enter irregular breathing
                
                % change state of most recent line
                fmS(1) = 4;
                
                % the previous state might have also bee IRR, check to see and
                % fix if so
                if fmS(1) == fmS(2) && fmIndStart(2) > 1
                    
                    % update beginning times, indices, statesof the most
                    % recent lines, starting with the most recent one
                    % also knots
                    fmTStart    = circshift(fmTStart,-1);
                    fmIndStart  = circshift(fmIndStart,-1);
                    fmS         = circshift(fmS,-1);
                    xKnot       = circshift(xKnot,-1);
                    
                end
                
            end
            
        end
        
    end
    
    % only do the following checks if we are not in the very first line
    if fmIndStart(1) > 1
        
        % check if current starting index is same as next most
        % recent (to within 133 ms)
        if fmTStart(1)-fmTStart(2) <= 2*cSLength
            %% local pruning
            
            [fmTStart,fmIndStart,fmS,xKnot] = matRad_FSMlocalPruning(fmTStart,fmIndStart,fmS,xKnot,velocityRangeS,dV);
            
        end
        
    end
    
end

[fmTStart,fmIndStart,xKnot,state] = matRad_FSMfinalization(x,t,deltaT,fmTStart,fmIndStart,fmS,xKnot,state);


end




