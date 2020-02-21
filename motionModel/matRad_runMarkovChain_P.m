function [l_simulated,relProb]= matRad_runMarkovChain_P(Pij,nSteps,l_init,maxProb)

if nargin < 4
    maxProb = false;
end

l_simulated = zeros(nSteps,1);

% the state of the first step in the sequence is given by l_init
step = 1;
l_simulated(step) = l_init;

if nargout > 1
    relProb = 1;
end

for step = 2:nSteps
    
    % extra current phase
    l_curr = l_simulated(step-1);
    
    if maxProb || nargout > 1
        % next phase is the maximum probability
        
        % maximum probability
        l_nextMax = find(Pij(l_curr,:) == max(Pij(l_curr,:)));
        l_nextMax = l_nextMax(randi(numel(l_nextMax)));
        
        % next phase is the maximum probability
        l_next = l_nextMax;
    else
        % do normal Monte Carlo
        
        % calculate cumulative probability distribution
        cumProb = cumsum(Pij(l_curr,:));
        
        % generate random number between 0 and 1
        r = rand;
        
        % determine next phase
        % Monte Carlo
        l_next = find(r < cumProb,1,'first');
    end
    l_simulated(step) = l_next;
    
    if nargout > 1
        relProb = relProb.*Pij(l_curr,l_next)./Pij(l_curr,l_nextMax);
    end
    
end

end

