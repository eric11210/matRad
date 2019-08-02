function motionModel = matRad_prepModelForOpt(prop4D)

% check if motion model already exists
if isfield(prop4D,'motionModel')
    
    % use model if it exists
    motionModel        = prop4D.motionModel;
    motionModel.type   = 'Markov';
    
    % strip model of the "out of bounds" phases
    motionModel = matRad_stripMarkovOOB(motionModel);
    
    % diagonalize matrix
    [motionModel.qij_V,motionModel.qij_D] = eig(motionModel.qij);
    
    % determine initial position phase
    posPhaseProb = accumarray(motionModel.indices.subPhase2PosPhase,motionModel.Pi_deltaTSample);
    initPosPhase = find(posPhaseProb == max(posPhaseProb));
    
    % determine initial probabililty
    motionModel.initProb = zeros(1,motionModel.indices.nSubPhases);
    % trigger on first phase and EOE
    initSubPhases   = motionModel.indices.subPhase2PosPhase == initPosPhase & motionModel.indices.subPhase2FS == 3;
    % let all subphases corresponding to the initial position phase
    % and FS (EOE) have the same probability
    motionModel.initProb(initSubPhases) = 1;
    % let all other subphases have the same (much smaller, but
    % nonzero) probability
    motionModel.initProb(~initSubPhases) = 1e-8;
    % normalize
    motionModel.initProb = motionModel.initProb./sum(motionModel.initProb);
    
else
    
    % if it doesn't, throw an error
    error('No motion model found in structure pln.propOpt.prop4d');
end

end

