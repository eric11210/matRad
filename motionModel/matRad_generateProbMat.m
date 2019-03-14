function model = matRad_generateProbMat(data)


% get variables from data struct
l_sample        = data.l_sample;
nSubPhases      = data.indices.nSubPhases;
deltaT_sample   = data.deltaT_sample;

% count number of transitions from one phase to the next
zij_deltaTSample = zeros(nSubPhases,nSubPhases);

for i = 1:numel(l_sample)
    
    if i == 1
        currPhase = l_sample(i);
        nextPhase = l_sample(i+1);
    elseif i == numel(l_sample)
        %prevPhase = p_sample(i-1);
        currPhase = l_sample(i);
    else
        %prevPhase = p_sample(i-1);
        currPhase = l_sample(i);
        nextPhase = l_sample(i+1);
    end
    
    if i ~= numel(l_sample)
        % if this isn't the last point, increment count by 1
        zij_deltaTSample(currPhase,nextPhase) = zij_deltaTSample(currPhase,nextPhase)+1;
    else
        % if it is, check if the current phase goes anywhere
        if ~any(zij_deltaTSample(currPhase,:))
            % if not, then let it go to the last phase with probability 1
            prevPhase = l_sample(i-1);
            
            zij_deltaTSample(currPhase,prevPhase) = zij_deltaTSample(currPhase,prevPhase)+1;
        end
    end
end

% construction of probability transition matrix
zi_deltaTSample = repmat(sum(zij_deltaTSample,2),1,nSubPhases);
Pij_deltaTSample = zij_deltaTSample./zi_deltaTSample;
Pij_deltaTSample(zi_deltaTSample == 0) = 0;

% construction of probability vector (time-homogeneous)
Pi_deltaTSample = zeros(nSubPhases,1);

for i = 1:nSubPhases
    Pi_deltaTSample(i) = nnz(l_sample == i)./numel(l_sample);
end

% construct transition rate matrix
qij = (Pij_deltaTSample-eye(nSubPhases))./deltaT_sample;

if false
    % determine possible transitions
    transij = zeros(nPhases,nPhases);
    for iSubPhase = 1:nSubPhases
        
        iPhase = subPhase2Phase(iSubPhase);
        if iPhase > nPhases
            continue
        end
        
        transij_inclExt = accumarray(subPhase2Phase,Pij_deltaTSample(iSubPhase,:));
        transij(iPhase,:) = transij(iPhase,:)+transij_inclExt(1:nPhases)';
    end
    transij = transij./repmat(nSubPhasePerPhase(1:nPhases),1,nPhases);
    
    transij(transij < 0.01) = 0;
    transij = logical(transij);
    % NOTE THAT THIS TRANSIJ IS FOR deltaT_sample
    % SHOULD DO MATRIX MULT TO GET A BETTER IDEA OF ACTUAL POSSIBLE TRANSITIONS
    % FROM ONE ANGLE TO THE NEXT
    % MAYBE NEED TO MOVE THIS TO THE PLN STRUCT
end

% put variables in model struct
model.Pij_deltaTSample  = Pij_deltaTSample;
model.Pi_deltaTSample   = Pi_deltaTSample;
model.qij               = qij;
model.deltaT_sample     = data.deltaT_sample;
model.indices           = data.indices;

end

