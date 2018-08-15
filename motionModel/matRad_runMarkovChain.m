function [l_simulated,relProb]= matRad_runMarkovChain(Pij,nSteps,l_init,maxProb)

if nargin < 4
    maxProb = false;
end

l_simulated = zeros(nSteps,1);

% the state of the first step in the sequence is given by l_init
step = 1;
l_simulated(step) = l_init;

relProb = 1;

for step = 2:nSteps
    
    % extra current phase
    l_curr = l_simulated(step-1);
    
    % calculate cumulative probability distribution
    cumProb = cumsum(Pij(l_curr,:));
    
    % generate random number between 0 and 1
    r = rand;
    
    % determine next phase
    % Monte Carlo
    l_nextRand = find(r < cumProb,1,'first');
    % maximum probability
    l_nextMax = find(Pij(l_curr,:) == max(Pij(l_curr,:)));
    l_nextMax = l_nextMax(randi(numel(l_nextMax)));
    
    if maxProb
        % next phase is the maximum probability
        l_next = l_nextMax;
    else
        % do normal Monte Carlo
        l_next = l_nextRand;
    end
    l_simulated(step) = l_next;
    
    relProb = relProb.*Pij(l_curr,l_next)./Pij(l_curr,l_nextMax);
    
end

end

