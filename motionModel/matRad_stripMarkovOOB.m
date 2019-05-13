function model = matRad_stripMarkovOOB(model)

% find the subphases which correspond to out of bounds position phases
indOOB_subPhase     = model.indices.subPhase2PosPhase == model.indices.nPosPhases-1 | model.indices.subPhase2PosPhase == model.indices.nPosPhases;
indOOB_posSubPhase  = model.indices.posSubPhase2PosPhase == model.indices.nPosPhases-1 | model.indices.posSubPhase2PosPhase == model.indices.nPosPhases;

% get new number of subPhases
nSubPhases_new = nnz(~indOOB_subPhase);

% delete indices in Pij, qij, Pi
model.Pij_deltaTSample(:,indOOB_subPhase)   = [];
model.Pij_deltaTSample(indOOB_subPhase,:)   = [];
model.qij(:,indOOB_subPhase)                = [];
model.qij(indOOB_subPhase,:)                = [];
model.Pi_deltaTSample(indOOB_subPhase)      = [];

% update indices (by deleting some)
model.indices.subPhase2PosSubPhase(indOOB_subPhase) = [];
model.indices.subPhase2VelSubPhase(indOOB_subPhase) = [];
model.indices.subPhase2State(indOOB_subPhase)       = [];
model.indices.subPhase2PosPhase(indOOB_subPhase)    = [];
model.indices.subPhase2VelPhase(indOOB_subPhase)    = [];
model.indices.subPhase2Phase(indOOB_subPhase)       = [];

model.indices.nSubPhasePerPhase     = accumarray(model.indices.subPhase2Phase,1);
model.indices.nSubPhasePerVelPhase  = accumarray(model.indices.subPhase2VelPhase,1);
model.indices.nSubPhasePerPosPhase  = accumarray(model.indices.subPhase2PosPhase,1);
model.indices.nSubPhases            = nSubPhases_new;

model.indices.posSubPhase2Pos(indOOB_posSubPhase)       = [];
model.indices.posSubPhase2PosPhase(indOOB_posSubPhase)  = [];
model.indices.nPosSubPhases                             = numel(model.indices.posSubPhase2PosPhase);
model.indices.nPosPhases                                = model.indices.nPosPhases-2;
model.indices.posPhase2Pos((end-1):end)                 = [];

%% continue cleaning the matrices until we don't have any more columns which sum to 0

% find columns summing to 0
normPij     = sum(model.Pij_deltaTSample,2);
deleteSubPhase   = normPij == 0;

while any(deleteSubPhase)
    
    % insert deleteInd in model struct
    model.deleteSubPhase = deleteSubPhase;
    
    % clean up matrices
    model = matRad_cleanProbMat(model);
    
    % find columns summing to 0
    normPij     = sum(model.Pij_deltaTSample,2);
    deleteSubPhase   = normPij == 0;
end

% renormalize Pij, qij, Pi
model.Pij_deltaTSample  = model.Pij_deltaTSample./normPij;
model.Pi_deltaTSample   = model.Pi_deltaTSample./sum(model.Pi_deltaTSample);
model.qij               = (model.Pij_deltaTSample-eye(size(model.Pij_deltaTSample)))./model.deltaT_sample;

end