function model = matRad_generateProbMat(data)


% get variables from data struct
l_sample        = data.l_sample;
nSubPhases      = data.indices.nSubPhases;
deltaT_sample   = data.deltaT_sample;
filtFactor      = data.filtFactor;

% count number of transitions from one phase to the next
zij_deltaTSample = zeros(nSubPhases,nSubPhases);

for i = 1:numel(l_sample)
    
    if i <= numel(l_sample)-filtFactor
        % normal current and next phase
        currPhase = l_sample(i);
        nextPhase = l_sample(i+filtFactor);
    else
        % let the next phase actually be the previous phase to prevent
        % normalization issues
        currPhase = l_sample(i);
        nextPhase = l_sample(i-filtFactor);
    end
    
    % if this isn't the last point, increment count by 1
    zij_deltaTSample(currPhase,nextPhase) = zij_deltaTSample(currPhase,nextPhase)+1;
end

% construction of probability transition matrix
zi_deltaTSample = repmat(sum(zij_deltaTSample,2),1,nSubPhases);
Pij_deltaTSample = zij_deltaTSample./zi_deltaTSample;
Pij_deltaTSample(zi_deltaTSample == 0) = 0;

% construction of probability vector (time-homogeneous)
Pi_deltaTSample = zeros(nSubPhases,1);

% prepare to delete subphases which are not reached
deleteSubPhase = zeros(nSubPhases,1,'logical');

for i = 1:nSubPhases
    Pi_deltaTSample(i) = nnz(l_sample(1:(end-filtFactor)) == i)./numel(l_sample(1:(end-filtFactor)));
    
    deleteSubPhase(i) = sum(Pij_deltaTSample(i,:)) == 0;
end

% construct transition rate matrix
identity = eye(nSubPhases);
identity(deleteSubPhase,deleteSubPhase) = 0;
qij = (Pij_deltaTSample-identity)./(deltaT_sample);

% put variables in model struct
model.Pij_deltaTSample  = Pij_deltaTSample;
model.Pi_deltaTSample   = Pi_deltaTSample;
model.qij               = qij;
model.deltaT_sample     = data.deltaT_sample;
model.indices           = data.indices;
model.deleteSubPhase    = deleteSubPhase;

end

