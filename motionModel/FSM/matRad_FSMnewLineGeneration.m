function [fmTStart,fmIndStart,xKnot,fmS] = matRad_FSMnewLineGeneration(x,t,deltaT,fmTStart,fmIndStart,fmS,xKnot,i,K,fV,dV,fSNext)

% update beginning times of the most recent lines,
% starting with the least recent one
fmTStart    = circshift(fmTStart,1);
fmTStart(1) = matRad_FSMfreeSplineFit(x(fmIndStart(1):i),t(fmIndStart(1):i),K,fV,dV);

% update knot positions
xKnot = circshift(xKnot,1);

% update starting indices
fmIndStart = ceil(fmTStart./deltaT);

% update states of the most recent lines, starting with
% the least recent one
fmS     = circshift(fmS,1);
fmS(1)  = fSNext;

end

