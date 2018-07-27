function tmm = generateTmm(x,t,nPhases,nSubPerPhase)

% convert raw data to phase
[l_sample,t_sample,deltaT_sample,subPhase2Phase,nSubPhases] = motionData2Phase(x,t,nPhases,nSubPerPhase);
nSubPhasePerPhase = accumarray(subPhase2Phase,1);

% construct probability matrix
[Pij_deltaTSample, Pi_deltaTSample] = generateProbMatrix(l_sample,nSubPhases);

%% move next two to generateProbMatrix?

% construct transition rate matrix
qij = (Pij_deltaTSample-eye(nSubPhases))./deltaT_sample;

% determine possible transitions
if true
    
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
else
    
    transij = zeros(nPhases+2,nPhases+2);
    for iSubPhase = 1:nSubPhases
        
        iPhase = subPhase2Phase(iSubPhase);
        if iPhase > nPhases
            %continue
        end
        
        transij_inclExt = accumarray(subPhase2Phase,Pij_deltaTSample(iSubPhase,:));
        transij(iPhase,:) = transij(iPhase,:)+transij_inclExt';
    end
    transij = transij./repmat(nSubPhasePerPhase,1,nPhases+2);
end

transij(transij < 0.01) = 0;
transij = logical(transij);
% NOTE THAT THIS TRANSIJ IS FOR deltaT_sample
% SHOULD DO MATRIX MULT TO GET A BETTER IDEA OF ACTUAL POSSIBLE TRANSITIONS
% FROM ONE ANGLE TO THE NEXT
% MAYBE NEED TO MOVE THIS TO THE PLN STRUCT


% generate tmm struct
tmm.qij                 = qij;
tmm.Pij_deltaTSample    = Pij_deltaTSample;
tmm.Pi_deltaTSample     = Pi_deltaTSample;
tmm.deltaT_sample       = deltaT_sample;
tmm.subPhase2Phase      = subPhase2Phase;
tmm.nSubPhasePerPhase   = nSubPhasePerPhase;
tmm.nSubPhases          = nSubPhases;
tmm.transij             = transij;

end