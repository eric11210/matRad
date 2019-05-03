function [lSimulated,tSimulated]= matRad_runMarkovChain_Q(qij,initProb,tStop)
% simulate a history of the Markov Chain

numFrames = numel(initProb);

% get rate parameters
lambda_vec = -diag(qij);

% get mean parameters
mu_vec = 1./lambda_vec;

% get probability transition matrix
% CHECK THIS
pij = qij.*mu_vec;
% take transpose to make memory access faster
pij_transpose = pij';

% estimate maximum number of transitions in time tStop
maxTrans = 10.*ceil(tStop.*lambda_vec);

% allocate frames and transition times
lSimulated = zeros(maxTrans,1);
tTrans     = zeros(maxTrans,1);

% sample first frame
lSimulated(1) = randsample(numFrames,1,true,initProb);

% initialize total time
tTot = 0;

% initialize step counter
step = 1;

% stop loop once total simulated exceeds stopping time
while tTot < tStop
    
    % get current frame
    lCurrent = lSimulated(step);
    
    % get exponential rate parameter
    mu = mu_vec(lCurrent);
    
    % sample from exponential distribution with mu parameter
    tSample = exprnd(mu);
    
    % insert sampled time into transition times
    tTrans(step) = tSample;
    
    % increase total time
    tTot = tTot+tSample;
    
    % sample next frame
    lSimulated(step+1) = randsample(numFrames,1,true,pij_transpose(:,lCurrent));
    
    % increase step
    step = step+1;
    
end

% delete any elements not used in lSimulated and tTrans
lSimulated((step+1):end)    = [];
tTrans(step:end)            = [];

% accumulate transition times to get simulated times
tSimulated = [0; cumsum(tTrans)];

end

