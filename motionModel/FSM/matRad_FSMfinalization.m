function [fmTStart,fmIndStart,xKnot,state] = matRad_FSMfinalization(x,t,deltaT,fmTStart,fmIndStart,fmS,xKnot,state)

% determine total number of spline fits to do
numSplineFits = find(xKnot ~= 0,1,'first')-1;

for splineFit = 1:numSplineFits
    
    % determine index of last knot position and time
    indLast = numSplineFits-splineFit+2;
    
    % determine last knot position and time
    xKnotLast = xKnot(indLast);
    tKnotLast = fmTStart(indLast);
    
    % do a partial spline fit
    if splineFit == numSplineFits
        [fmTStartNew,xKnotNew] = matRad_FSMpartialSplineFit(x(fmIndStart(2):(fmIndStart(1)-1)),x(fmIndStart(1):end),t(fmIndStart(2):(fmIndStart(1)-1)),t(fmIndStart(1):end),xKnotLast,tKnotLast);
    else
        [fmTStartNew,xKnotNew] = matRad_FSMpartialSplineFit(x(fmIndStart(indLast):(fmIndStart(indLast-1)-1)),x(fmIndStart(indLast-1):(fmIndStart(indLast-2)-1)),t(fmIndStart(indLast):(fmIndStart(indLast-1)-1)),t(fmIndStart(indLast-1):(fmIndStart(indLast-2)-1)),xKnotLast,tKnotLast);
    end
    
    % update starting times
    fmTStart(indLast-1) = fmTStartNew;
    
    % update knot positions
    xKnot(indLast-1) = xKnotNew;
    
    % update starting indices
    fmIndStart = ceil(fmTStart./deltaT);
    
    % enter finalized state
    state(fmIndStart(indLast):(fmIndStart(indLast-1)-1)) = fmS(indLast);
    
end

% enter finalized states for the last line
state(fmIndStart(1):end) = fmS(1);

end

