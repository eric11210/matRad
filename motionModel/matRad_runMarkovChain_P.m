function lSimulated = matRad_runMarkovChain_P(motionModel,nSteps,maxProb)
% simulate a history of the Markov Chain using the P matrix

% extract transition matrix and initial probabilities
Pij             = motionModel.Pij_deltaTSample;
initProb        = motionModel.initProb;

if nargin < 3
    maxProb = false;
end

% allocate states
lSimulated = zeros(nSteps,1);

% determine number of states in Markov chain
numStates = numel(initProb);

% sample first frame
lSimulated(1) = randsample(numStates,1,true,initProb);

% loop through the steps
for step = 2:nSteps
    
    % get current state
    lCurrent = lSimulated(step-1);
    
    if maxProb
        % next phase is the maximum probability
        
        % maximum probability
        lNext = find(Pij(lCurrent,:) == max(Pij(lCurrent,:)));
        lNext = lNext(randi(numel(lNext)));
        
    else
        % do normal Monte Carlo
        
        % sample from current row of Pij
        lNext = randsample(numStates,1,true,Pij(lCurrent,:));
    end
    
    % insert next state
    lSimulated(step) = lNext;
end

end

