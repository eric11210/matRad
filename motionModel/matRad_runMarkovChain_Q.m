function [lSimulated,tSimulated] = matRad_runMarkovChain_Q(motionModel,tTrans,maxProb)
% simulate a history of the Markov Chain using the Q matrix

% extract transition matrix and initial probabilities
%qij         = motionModel.qij;
qij_D       = motionModel.qij_D;
qij_V       = motionModel.qij_V;
initProb    = motionModel.initProb;

if nargin < 3
    maxProb = false;
end

% determine number of states in Markov chain
numStates = numel(initProb);

% determine number of steps to take
nSteps = numel(tTrans);

% allocate frames
lSimulated = zeros(nSteps+1,1);

% accumulate transition times to get simulated times
tSimulated = [0; cumsum(tTrans)];

% sample first frame
lSimulated(1) = randsample(numStates,1,true,initProb);

% loop through the steps
for step = 1:nSteps
    
    % get current state
    lCurrent = lSimulated(step);
    
    % get time since last step
    tStep = tTrans(step);
    
    % calculate transition probability matrix for transition time tTrans
    %Pij_tStep = expm(qij.*tStep);
    Pij_tStep      = real(qij_V*diag(exp(diag(qij_D.*tStep)))/qij_V);
    
    % set any negative values to 0 and renormalize
    Pij_tStep(Pij_tStep < 0)    = 0;
    Pij_tStep                   = Pij_tStep./sum(Pij_tStep,2);
    
    % sample next state
    if maxProb
        [~,lNext]   = max(Pij_tStep(lCurrent,:));
    else
        lNext = randsample(numStates_real,1,true,Pij_tStep(lCurrent,:));
    end
    
    % insert next state
    lSimulated(step+1) = lNext;
    
end

end

