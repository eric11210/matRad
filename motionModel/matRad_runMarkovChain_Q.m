function [lSimulated,tSimulated]= matRad_runMarkovChain_Q(qij,initProb,tTrans,indexMap)
% simulate a history of the Markov Chain

if nargin < 4
    indexMap = (1:numel(initProb))';
end

norm = accumarray(indexMap,1);

[indexMap_gridJ,indexMap_gridI] = meshgrid(indexMap);

numFrames = numel(initProb);

numPhases = max(indexMap);

% determine number of steps to take
nSteps = numel(tTrans);

% allocate frames
lSimulated = zeros(nSteps+1,1);

% accumulate transition times to get simulated times
tSimulated = [0; cumsum(tTrans)];

% sample first frame
lSimulated(1) = indexMap(randsample(numFrames,1,true,initProb));

% loop through the steps
for step = 1:nSteps
    
    %% simulation
    
    % get current frame
    lCurrent = lSimulated(step);
    
    % get time since last step
    tStep = tTrans(step);
    
    % calculate transition probability to arrive at i at T and derivative
    Pij_tStep = expm(qij.*tStep);
    
    % map the Pij matrix
    PijMapped_tStep = accumarray([indexMap_gridI(:) indexMap_gridJ(:)],Pij_tStep(:))./norm;
    
    % sample next frame
    lNext = indexMap(randsample(numPhases,1,true,PijMapped_tStep(lCurrent,:)));
    
    % insert next frame
    lSimulated(step+1) = lNext;
    
end

end

