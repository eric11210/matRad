function [Pij_transT,Pij_transT_dot,Pi_T,Pi_T_dot] = matRad_transAndTProb_GPU(transT,T,motionModel)

persistent qij_GPU initProb_GPU Pij_transT_GPU Pij_transT_dot_GPU Pi_T_GPU Pi_T_dot_GPU

% calculate probabilities differently depending on type of motion model
switch motionModel.type
    case 'Markov'
        
        % extract transition matrix and initial probabilities
        qij         = motionModel.qij;
        initProb    = motionModel.initProb;
        
        if isempty(qij_GPU)
            qij_GPU = gpuArray(qij);
            initProb_GPU = gpuArray(initProb);
        end
        
        % calculate transition probability matrix and derivative
        Pij_transT_GPU       = expm(qij_GPU.*transT);
        Pij_transT_dot_GPU   = qij_GPU*Pij_transT_GPU;
        
        % calculate probability to arrive at i at T and derivative
        Pi_T_GPU     = initProb_GPU*expm(qij_GPU.*T);
        Pi_T_dot_GPU = Pi_T_GPU*qij_GPU;
        
        % gather
        Pij_transT = gather(Pij_transT_GPU);
        Pij_transT_dot = gather(Pij_transT_dot_GPU);
        Pi_T = gather(Pi_T_GPU);
        Pi_T_dot = gather(Pi_T_dot_GPU);
        
        
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