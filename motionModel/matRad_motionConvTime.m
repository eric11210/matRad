function model = matRad_motionConvTime(model,options)

% extract variables from struct
Pij_deltaTSample    = model.Pij_deltaTSample;
Pi_deltaTSample     = model.Pi_deltaTSample;
deltaT_sample       = model.deltaT_sample;
subPhase2Phase      = model.indices.subPhase2Phase;
nSubPhases          = model.indices.nSubPhases;
nSubPhasePerPhase   = model.indices.nSubPhasePerPhase;
percRMSD_targ       = options.percRMSD_targ;

Pi_recombined = accumarray(subPhase2Phase,Pi_deltaTSample);

% initial probability vector
initPVec = zeros(1,nSubPhases);
initPVec(subPhase2Phase == 1) = 1./nSubPhasePerPhase(1);

% initial error
convergeT = 0;
Pij_T = eye(nSubPhases,nSubPhases);
Pi_T_recombined = accumarray(subPhase2Phase,initPVec*Pij_T);
SD = ((Pi_T_recombined-Pi_recombined)./Pi_recombined).^2;
SD(Pi_recombined == 0) = 0;
percRMSD = 100*sqrt(mean(SD));

while percRMSD > percRMSD_targ
    
    % take another step
    convergeT = convergeT+deltaT_sample;
    
    % matrix mult
    Pij_T = Pij_T*Pij_deltaTSample;
    
    % sum over subphases
    Pi_T_recombined = accumarray(subPhase2Phase,initPVec*Pij_T);
    
    % calculate error
    SD = ((Pi_T_recombined-Pi_recombined)./Pi_recombined).^2;
    SD(Pi_recombined == 0) = 0;
    percRMSD = 100*sqrt(mean(SD));
    
end

% put variables in struct
model.convergeT = convergeT;

end

