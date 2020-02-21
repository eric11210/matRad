function model = matRad_cleanProbMat(model)

% get subPhases to be deleted, remove subfield
deleteSubPhase  = model.deleteSubPhase;
model           = rmfield(model,'deleteSubPhase');

% delete indices in Pij, qij, Pi
model.Pij_deltaTSample(:,deleteSubPhase)    = [];
model.Pij_deltaTSample(deleteSubPhase,:)    = [];
model.qij(:,deleteSubPhase)                 = [];
model.qij(deleteSubPhase,:)                 = [];
model.Pi_deltaTSample(deleteSubPhase)       = [];

% diagonalize matrix
[model.qij_V,model.qij_D]   = eig(model.qij);

% get new number of subPhases
nSubPhases_new = nnz(~deleteSubPhase);

% update indices (by deleting some)
model.indices.subPhase2PosSubPhase(deleteSubPhase)  = [];
model.indices.subPhase2VelSubPhase(deleteSubPhase)  = [];
model.indices.subPhase2State(deleteSubPhase)        = [];
model.indices.subPhase2FS(deleteSubPhase)           = [];
model.indices.subPhase2PosPhase(deleteSubPhase)     = [];
model.indices.subPhase2VelPhase(deleteSubPhase)     = [];
model.indices.subPhase2Phase(deleteSubPhase)        = [];

model.indices.nSubPhasePerPhase     = accumarray(model.indices.subPhase2Phase,1);
model.indices.nSubPhasePerVelPhase  = accumarray(model.indices.subPhase2VelPhase,1);
model.indices.nSubPhasePerPosPhase  = accumarray(model.indices.subPhase2PosPhase,1);
model.indices.nSubPhases            = nSubPhases_new;

end