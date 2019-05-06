function [Pij_transT,Pij_transT_dot,Pi_T,Pi_T_dot] = matRad_transAndTProb(transT,T,motionModel)

% calculate probabilities differently depending on type of motion model
switch motionModel.type
    case 'Markov'
        
        % extract transition matrix and initial probabilities
        qij         = motionModel.qij;
        initProb    = motionModel.initProb;
        
        % calculate transition probability matrix and derivative
        Pij_transT      = expm(qij.*transT);
        Pij_transT_dot  = qij*Pij_transT;
        
        % calculate probability to arrive at i at T and derivative
        Pi_T        = initProb*expm(qij.*T);
        Pi_T_dot    = Pi_T*qij;
        
        % [subPhase2PosPhase_gridJ, subPhase2PosPhase_gridI] = meshgrid(subPhase2PosPhase);
        
        % norm = repmat(accumarray(subPhase2PosPhase_gridI(:),Pij(:)),[1 nPosPhases];
        
        % PIJ = accumarray([subPhase2PosPhase_gridI(:) subPhase2PosPhase_gridJ(:)],Pij(:))./norm;
        
        % PIJ(isnan(PIJ)) = 0;
        
        % PI = accumarray(subPhase2PosPhase,Pi);
        
    case 'single'
        
        % allocate the probability matrices and vectors
        Pij_transT  = zeros(motionModel.numPhases,motionModel.numPhases);
        Pi_T        = zeros(1,motionModel.numPhases);
        
        % determine phases at beginning and end of arc
        phase_I = motionModel.lSimulated(motionModel.tSimulated <= T);
        phase_F = motionModel.lSimulated(motionModel.tSimulated <= (T+transT));
        phase_I = phase_I(end);
        phase_F = phase_F(end);
        
        % use these phases to construct the probability matrices and vectors
        Pij_transT(phase_I,phase_F) = 1;
        Pi_T(phase_I)               = 1;
        
        % the derivatives can just be equal to 0
        Pij_transT_dot  = zeros(motionModel.numPhases,motionModel.numPhases);
        Pi_T_dot        = zeros(1,motionModel.numPhases);
        
end

end