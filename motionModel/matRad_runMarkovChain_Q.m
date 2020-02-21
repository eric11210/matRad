function [lSimulated,tSimulated] = matRad_runMarkovChain_Q(motionModel,tTrans,indexMap,maxProb)
% simulate a history of the Markov Chain

% extract transition matrix and initial probabilities
%qij         = motionModel.qij;
qij_D       = motionModel.qij_D;
qij_V       = motionModel.qij_V;
initProb    = motionModel.initProb;

if nargin < 3
    indexMap = (1:numel(initProb))';
    maxProb = false;
end

if nargin < 4
    maxProb = false;
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
    
    % calculate transition probability to arrive at i at T
    %Pij_tStep = expm(qij.*tStep);
    Pij_tStep      = real(qij_V*diag(exp(diag(qij_D.*tStep)))/qij_V);
    
    % set any negative values to 0 and renormalize
    Pij_tStep(Pij_tStep < 0)    = 0;
    Pij_tStep                   = Pij_tStep./sum(Pij_tStep,2);
    
    % map the Pij matrix
    PijMapped_tStep = accumarray([indexMap_gridI(:) indexMap_gridJ(:)],Pij_tStep(:))./norm;
    
    % sample next frame
    if maxProb
        [~,lNext]   = max(PijMapped_tStep(lCurrent,:));
        lNext       = indexMap(lNext);
    else
        lNext = indexMap(randsample(numPhases,1,true,PijMapped_tStep(lCurrent,:)));
    end
    
    % insert next frame
    lSimulated(step+1) = lNext;
    
end

end

