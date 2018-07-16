function convergeT = convergenceTime(Pij_deltaTSample,Pi_deltaTSample,initPVec,subsForAccum,targetRelError,deltaT_sample,tInit)

Pi_recombined = accumarray(subsForAccum,Pi_deltaTSample);

if nargin < 7
    nInit = 0;
else
    nInit = max(floor(tInit./deltaT_sample),0);
end

% initial error
nStep = nInit;
Pn_recombined = accumarray(subsForAccum,initPVec*Pij_deltaTSample^nStep);
relError = sum(((Pn_recombined-Pi_recombined)./Pi_recombined).^2);

while relError > targetRelError
    
    % take another step
    nStep = nStep+1;
    % recalculate error
    Pn_recombined = accumarray(subsForAccum,initPVec*Pij_deltaTSample^nStep);
    relError = sum(((Pn_recombined-Pi_recombined)./Pi_recombined).^2);
end

convergeT = nStep*deltaT_sample;

end