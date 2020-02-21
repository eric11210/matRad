function model = matRad_motionConvTime(model,options)

%% setup

% extract model variables from struct
Pij_deltaTSample        = model.Pij_deltaTSample;
deltaT_sample           = model.deltaT_sample;
subPhase2PosPhase       = model.indices.subPhase2PosPhase;
nPosPhases              = model.indices.nPosPhases;
nSubPhases              = model.indices.nSubPhases;
nSubPhasePerPosPhase    = model.indices.nSubPhasePerPhase;
l2_targ                 = options.l2_targ;

%% determine steady state of distribution & decay timescale

% get eigenvalues and eigenvectors of Pij
[r,d,l] = eig(Pij_deltaTSample);

% calculate steady state of the phase probability distribution
idEVec_ind          = find(abs(diag(d)-1) < 1e-10);
Pi_inf              = l(:,idEVec_ind)./sum(l(:,idEVec_ind));
Pi_inf_recombined   = accumarray(subPhase2PosPhase,Pi_inf,[nPosPhases 1]);

% get decay timescale, stepsize
dAbsVec                         = sort(diag(abs(d)),'descend');
dAbsVec(abs(dAbsVec-1) < 1e-10) = [];
decayT                          = round(-1./log(dAbsVec(2)));
stepSize                        = round(decayT/8);

%% estimate convergence time

% pre-calculate matrix power
Pij_deltaTSample_pStepsize = Pij_deltaTSample^stepSize;

% initial probability vector
initPVec = zeros(1,nSubPhases);
initPVec(subPhase2PosPhase == 1) = 1./nSubPhasePerPosPhase(1);

% initial error
convergeT = 0;
Pij_T = eye(nSubPhases,nSubPhases);
Pi_T_recombined = accumarray(subPhase2PosPhase,initPVec*Pij_T,[nPosPhases 1]);
l2 = sqrt(sum((Pi_T_recombined-Pi_inf_recombined).^2));

while l2 > l2_targ && convergeT/decayT < 100
    
    % take another step
    convergeT = convergeT+stepSize;
    
    % matrix mult
    Pij_T = Pij_T*Pij_deltaTSample_pStepsize;
    
    % sum over subphases
    Pi_T_recombined = accumarray(subPhase2PosPhase,initPVec*Pij_T,[nPosPhases 1]);
    
    % calculate l2
    l2 = sqrt(sum((Pi_T_recombined-Pi_inf_recombined).^2));
    
end

% multiply deltaT into convergence time
model.convergeT = convergeT.*deltaT_sample;

end

