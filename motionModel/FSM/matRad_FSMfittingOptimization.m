function [fmTStart,fmIndStart,xKnot,state,cLambdaVec,cThetaVec,vEOEVec] = matRad_FSMfittingOptimization(x,t,deltaT,fmTStart,fmIndStart,fmS,xKnot,state,cLambdaVec,cThetaVec,vEOEVec)

% only do this if we have generated enough lines
if fmIndStart(4) > 1
    
    % first check if any of the starting indices are the same
    % if they are, bump the most recent one up by one
    if fmIndStart(4)-1 <= fmIndStart(5)
        fmIndStart(4) = fmIndStart(5)+2;
    end
    if fmIndStart(3)-1 <= fmIndStart(4)
        fmIndStart(3) = fmIndStart(4)+2;
    end
    if fmIndStart(2)-1 <= fmIndStart(3)
        fmIndStart(2) = fmIndStart(3)+2;
    end
    
    % check if this is our first knot that we are
    % generating
    if fmIndStart(5) == 1
        % first knot
        
        % do a spline fit (just don't pass any knot
        % argument to the function
        [fmTStartNew,xKnotNew] = matRad_FSMpartialSplineFit(x(fmIndStart(5):(fmIndStart(4)-1)),x(fmIndStart(4):(fmIndStart(3)-1)),t(fmIndStart(5):(fmIndStart(4)-1)),t(fmIndStart(4):(fmIndStart(3)-1)));
        
    else
        
        % determine last knot position and time
        xKnotLast = xKnot(5);
        tKnotLast = fmTStart(5);
        
        % do a partial spline fit
        [fmTStartNew,xKnotNew] = matRad_FSMpartialSplineFit(x(fmIndStart(5):(fmIndStart(4)-1)),x(fmIndStart(4):(fmIndStart(3)-1)),t(fmIndStart(5):(fmIndStart(4)-1)),t(fmIndStart(4):(fmIndStart(3)-1)),xKnotLast,tKnotLast);
        
    end
    
    % update starting times
    fmTStart(4) = fmTStartNew;
    
    % update knot positions
    xKnot(4) = xKnotNew;
    
    % update starting indices
    fmIndStart = ceil(fmTStart./deltaT);
    
    % enter finalized state
    state(fmIndStart(5):(fmIndStart(4)-1)) = fmS(5);
    
    % dynamic adaptation
    
    % check if the finalized state is either exhale or inhale
    if fmS(5) == 1 || fmS(5) == 2
        
        % update thresholds on amplitude, starting with least
        % recent
        cLambdaVec = circshift(cLambdaVec,1);
        % the most recent threshold is one third of the amplitude of the
        % finalized line
        cLambdaVec(1) = matRad_FSMgetAmplitude(x(fmIndStart(5):(fmIndStart(4)-1)),x(fmIndStart(5):(fmIndStart(4)-1)),deltaT)/3;
        
        % now check if the finalized state is exhale
        if fmS(5) == 2
            
            % update thresholds on velocity difference,
            % starting with least recent
            cThetaVec = circshift(cThetaVec,1);
            % the most recent threshold is one quarter of the
            % velocity of the finalized line
            cThetaVec(1) = matRad_FSMgetFitVelocity(x(fmIndStart(5):(fmIndStart(4)-1)),deltaT)./4;
            
        end
        
    elseif fmS(5) == 3
        
        vEOEVec = circshift(vEOEVec,1);
        vEOEVec(1) = matRad_FSMgetFitVelocity(x(fmIndStart(5):(fmIndStart(4)-1)),deltaT);
        
    end
    
end

end