function [lSimulated,tSimulated]= matRad_runMarkovChain_Q(qij,initProb,tTrans)
% simulate a history of the Markov Chain

numFrames = numel(initProb);

% determine number of steps to take
nSteps = numel(tTrans);

% allocate frames
lSimulated = zeros(nSteps+1,1);

% accumulate transition times to get simulated times
tSimulated = [0; cumsum(tTrans)];

% sample first frame
lSimulated(1) = randsample(numFrames,1,true,initProb);

% loop through the steps
for step = 1:nSteps
    
    %% simulation
    
    % get current frame
    lCurrent = lSimulated(step);
    
    % get time since last step
    tStep = tTrans(step);
    
    % calculate transition probability to arrive at i at T and derivative
    Pij_tStep = expm(qij.*tStep);
    
    % sample next frame
    lNext = randsample(numFrames,1,true,Pij_tStep(lCurrent,:));
    
    % insert next frame
    lSimulated(step+1) = lNext;
    
end

end

