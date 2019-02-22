function [fmTStart,fmIndStart,xKnot,fmS,state,cLambdaVec,cThetaVec,vEOEVec] = matRad_FSMeoeEnforcement(x,t,deltaT,fmTStart,fmIndStart,fmS,xKnot,state,cLambdaVec,cThetaVec,vEOEVec,i,K)

%% find the pattern of the most recent cycle
% start by finding indices when state is exhale
indInh = find(fmS == 2,2,'first');

%if indInh(1) == 1
    % if the most recent state is inhale, then the start of the cycle is the
    % second most recent inhale
%    indCycle = indInh(2);
%else
    % otherwise the start of the cycle is the most recent exhale
    indCycle = indInh(1);
%end

% determine pattern
patCycle = num2str(fmS(1:indCycle)');

%% do EOE enforcement

% what we do here depends on what the pattern is
switch patCycle
    
    case '1  3  2'
        
        % do a free spline fit from fmIndStart(1) to i
        % the assumption here is that part of this region is inhale, and
        % part is EOE; we want to find the break point
        
        % first define "d" to be half the distance between fmIndStart(1) and i
        d = round((i-fmIndStart(1)+1)./2);
        
        % now find average velocities of the two lines
        v1 = matRad_FSMgetFitVelocity(x(fmIndStart(1):(fmIndStart(1)+d-1)),deltaT);
        v2 = matRad_FSMgetFitVelocity(x((fmIndStart(1)+d):i),deltaT);
        
        % find start time of inhale
        fmTStart(1) = matRad_FSMfreeSplineFit(x(fmIndStart(1):i),t(fmIndStart(1):i),d,v1,v2);
        
        % update starting indices
        fmIndStart = ceil(fmTStart./deltaT);
        
        
        % now do a free spline fit from fmIndStart(3) to fmIndStart(1)-1
        % (the updated one)
        
        % first define "d" to be the distance between fmIndStart(2) and
        % fmIndStart(1)-1
        d = fmIndStart(1)-fmIndStart(2);
        
        % correct
        if fmIndStart(3)+d == fmIndStart(1)-1
            d = d-1;
        end
        
        % now find average velocities of the two lines
        v1 = matRad_FSMgetFitVelocity(x(fmIndStart(3):(fmIndStart(3)+d-1)),deltaT);
        v2 = matRad_FSMgetFitVelocity(x((fmIndStart(3)+d):(fmIndStart(1)-1)),deltaT);
        
        % find start time of EOE
        fmTStart(2) = matRad_FSMfreeSplineFit(x(fmIndStart(3):(fmIndStart(1)-1)),t(fmIndStart(3):(fmIndStart(1)-1)),d,v1,v2);
        
        % update starting indices
        fmIndStart = ceil(fmTStart./deltaT);
        
    case '1  4  2'
        
        % do a free spline fit from fmIndStart(1) to i
        % the assumption here is that part of this region is inhale, and
        % part is EOE; we want to find the break point
        
        % first define "K" to be half the distance between fmIndStart(1) and i
        d = round((i-fmIndStart(1)+1)./2);
        
        % now find average velocities of the two lines
        v1 = matRad_FSMgetFitVelocity(x(fmIndStart(1):(fmIndStart(1)+d-1)),deltaT);
        v2 = matRad_FSMgetFitVelocity(x((fmIndStart(1)+d):i),deltaT);
        
        % find start time of inhale
        fmTStart(1) = matRad_FSMfreeSplineFit(x(fmIndStart(1):i),t(fmIndStart(1):i),d,v1,v2);
        
        % update starting indices
        fmIndStart = ceil(fmTStart./deltaT);
        
        
        % now do a free spline fit from fmIndStart(2) to fmIndStart(1)-1
        % (the updated one)
        
        % first define "d" to be half the distance between fmIndStart(1) and i
        d = round((fmIndStart(1)-fmIndStart(2))./2);
        
        % now find average velocities of the two lines
        v1 = matRad_FSMgetFitVelocity(x(fmIndStart(2):(fmIndStart(2)+d-1)),deltaT);
        v2 = matRad_FSMgetFitVelocity(x((fmIndStart(2)+d):(fmIndStart(1)-1)),deltaT);
        
        % find start time of EOE
        fmTStart(2) = matRad_FSMfreeSplineFit(x(fmIndStart(2):(fmIndStart(1)-1)),t(fmIndStart(2):(fmIndStart(1)-1)),d,v1,v2);
        
        % update starting indices
        fmIndStart = ceil(fmTStart./deltaT);
        
        
        % now do a partial spline fit from the knot at fmIndStart(4) to
        % fmIndStart(2)-1
        
        % determine last knot position and time
        xKnotLast = xKnot(4);
        tKnotLast = fmTStart(4);
        
        % do a partial spline fit
        [fmTStartNew,xKnotNew] = matRad_FSMpartialSplineFit(x(fmIndStart(4):(fmIndStart(3)-1)),x(fmIndStart(3):(fmIndStart(2)-1)),t(fmIndStart(4):(fmIndStart(3)-1)),t(fmIndStart(3):(fmIndStart(2)-1)),xKnotLast,tKnotLast);
        
        % update starting times
        fmTStart(3) = fmTStartNew;
        
        % update knot positions
        xKnot(3) = xKnotNew;
        
        % update starting indices
        fmIndStart = ceil(fmTStart./deltaT);
        
        % since we've changed fmIndStart(3), we need to reenter finalized state
        state(fmIndStart(4):(fmIndStart(3)-1)) = fmS(4);
        
        % change fmS(2) to EOE
        fmS(2) = 3;
        
    case '2'
        
        % do a free spline fit from fmIndStart(1) to i-K
        % the assumption here is that part of this region is exhale, and
        % part is EOE; we want to find the break point
        
        % first define "d" to be half the distance between fmIndStart(1)
        % and i-K
        d = round((i-K-fmIndStart(1)+1)./2);
        
        % now find average velocities of the two lines
        v1 = matRad_FSMgetFitVelocity(x(fmIndStart(1):(fmIndStart+d-1)),deltaT);
        v2 = matRad_FSMgetFitVelocity(x((fmIndStart(1)+d):(i-K)),deltaT);
        
        % update times/knots/states
        fmTStart    = circshift(fmTStart,1);
        xKnot       = circshift(xKnot,1);
        fmS         = circshift(fmS,1);
        fmS(1)      = 3;
        
        % find start time of EOE
        fmTStart(1) = matRad_FSMfreeSplineFit(x(fmIndStart(1):(i-K)),t(fmIndStart(1):(i-K)),d,v1,v2);
        
        % update starting indices
        fmIndStart = ceil(fmTStart./deltaT);
        
        % fitting optimization
        [fmTStart,fmIndStart,xKnot,state,cLambdaVec,cThetaVec,vEOEVec] = matRad_FSMfittingOptimization(x,t,deltaT,fmTStart,fmIndStart,fmS,xKnot,state,cLambdaVec,cThetaVec,vEOEVec);
        
        
        % now do a free spline fit from fmIndStart(1) to i
        % (the updated one)
        
        % first define "d" to be half the distance between fmIndStart(1) and i
        d = round((i-fmIndStart(1)+1)./2);
        
        % now find average velocities of the two lines
        v1 = matRad_FSMgetFitVelocity(x(fmIndStart(1):(fmIndStart(1)+d-1)),deltaT);
        v2 = matRad_FSMgetFitVelocity(x((fmIndStart(1)+d):i),deltaT);
        
        % update times/knots/states
        fmTStart    = circshift(fmTStart,1);
        xKnot       = circshift(xKnot,1);
        fmS         = circshift(fmS,1);
        fmS(1)      = 1;
        
        % find start time of EOE
        fmTStart(1) = matRad_FSMfreeSplineFit(x(fmIndStart(1):i),t(fmIndStart(1):i),d,v1,v2);
        
        % update starting indices
        fmIndStart = ceil(fmTStart./deltaT);
        
        % fitting optimization
        [fmTStart,fmIndStart,xKnot,state,cLambdaVec,cThetaVec,vEOEVec] = matRad_FSMfittingOptimization(x,t,deltaT,fmTStart,fmIndStart,fmS,xKnot,state,cLambdaVec,cThetaVec,vEOEVec);
        
    case '1  4  3  2'
        
        error('nothing');
        
    otherwise
        
        % give error message when we haven't yet programmed behaviour for
        % this pattern
        error('Error when doing EOE enforcement: pattern %s not programmed.',patCycle);
        
end

end

