function [l_simulated,prob] = MCMC(Pij_sample,nSteps,l_init,maxProb)

if nargin < 4
    maxProb = false;
end

l_simulated = zeros(nSteps,1);

% the state of the first step in the sequence is given by l_init
step = 1;
l_simulated(step) = l_init;

prob = 1;

for step = 2:nSteps
    
    % extra current phase
    l_curr = l_simulated(step-1);
    
    % calculate cumulative probability distribution
    cumProb = cumsum(Pij_sample(l_curr,:));
    
    % generate random number between 0 and 1
    r = rand;
    
    % determine next phase
    if maxProb
        % next phase is the maximum probability
        l_next = find(Pij_sample(l_curr,:) == max(Pij_sample(l_curr,:)));
    else
        % do normal Monte Carlo
        l_next = find(r < cumProb,1,'first');
    end
    l_simulated(step) = l_next;
    
    prob = prob*Pij_sample(l_curr,l_next);
    
end

end