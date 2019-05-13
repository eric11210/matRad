function [Pij_transT,Pij_transT_dot,Pi_T,Pi_T_dot] = matRad_transAndTProb(transT,T,motionModel)

% calculate probabilities differently depending on type of motion model
switch motionModel.type
    case 'Markov'
        
        % extract transition matrix and initial probabilities
        qij         = motionModel.qij;
        initProb    = motionModel.initProb;
        
        % calculate transition probability matrix and derivative
        PijSub_transT       = expm(qij.*transT);
        PijSub_transT_dot   = qij*PijSub_transT;
        
        % calculate probability to arrive at i at T and derivative
        PiSub_T     = initProb*expm(qij.*T);
        PiSub_T_dot = PiSub_T*qij;
        
        % accumulate probabilities into position phases
        [subPhase2PosPhase_gridJ, subPhase2PosPhase_gridI] = meshgrid(motionModel.indices.subPhase2PosPhase);
        
        norm            = motionModel.indices.nSubPhasePerPosPhase;
        Pij_transT      = accumarray([subPhase2PosPhase_gridI(:) subPhase2PosPhase_gridJ(:)],PijSub_transT(:))./norm;
        Pij_transT_dot  = accumarray([subPhase2PosPhase_gridI(:) subPhase2PosPhase_gridJ(:)],PijSub_transT_dot(:))./norm;
        
        Pi_T        = accumarray(motionModel.indices.subPhase2PosPhase,PiSub_T);
        Pi_T_dot    = accumarray(motionModel.indices.subPhase2PosPhase,PiSub_T_dot);
        
    case 'single'
        
        % allocate the probability matrices and vectors
        Pij_transT  = zeros(motionModel.numPhases,motionModel.numPhases);
        Pi_T        = zeros(1,motionModel.numPhases);
        
        % determine phases at beginning and end of arc
        phase_I = motionModel.lSimulated(abs(motionModel.tSimulated-T) <= eps.*T);
        phase_F = motionModel.lSimulated(abs(motionModel.tSimulated-(T+transT)) <= eps.*(T+transT));
        
        if isempty(phase_I) || isempty(phase_F)
            error('Phases not found!');
        end
        
        % use these phases to construct the probability matrices and vectors
        Pij_transT(phase_I,phase_F) = 1;
        Pi_T(phase_I)               = 1;
        
        % the derivatives can just be equal to 0
        Pij_transT_dot  = zeros(motionModel.numPhases,motionModel.numPhases);
        Pi_T_dot        = zeros(1,motionModel.numPhases);
        
end

end