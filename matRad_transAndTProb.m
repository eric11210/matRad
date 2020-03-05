function [Pij_transT,Pij_transT_dot,Pi_T,Pi_T_dot] = matRad_transAndTProb(transT,T,motionModel)

% calculate probabilities differently depending on type of motion model
switch motionModel.type
    case 'Markov_Q'
        % Markov chain based on transition rate matrix (Q)
        % no longer used except for trivial, non-4D optimization
        
        % extract transition rate matrix and initial probabilities
        qij         = motionModel.qij;
        qij_D       = motionModel.qij_D;
        qij_V       = motionModel.qij_V;
        initProb    = motionModel.initProb;
        
        % calculate transition probability matrix and derivative
        %Pij_transT      = expm(qij.*transT);
        if transT == 0
            Pij_transT = eye(motionModel.indices.nSubPhases);
        else
            Pij_transT = abs(qij_V*diag(exp(diag(qij_D.*transT)))/qij_V);
        end
        Pij_transT_dot  = qij*Pij_transT;
        
        % calculate probability to arrive at i at T and derivative
        %Pi_T        = initProb*expm(qij.*T);
        if T == 0
            Pi_T = initProb;
        else
            Pi_T = initProb*abs(qij_V*diag(exp(diag(qij_D.*T)))/qij_V);
        end
        Pi_T_dot    = Pi_T*qij;
        
    case 'Markov_P'
        % Markov chain based on transition probability matrix (P)
        
        % extract transition matrix, deltaT, and initial probabilities
        Pij         = motionModel.Pij_deltaTSample;
        initProb    = motionModel.initProb;
        deltaT      = motionModel.deltaT_sample;
        
        % find nearest integer multiples of transT/deltaTSample and
        % T/deltaTSample
        transT_n    = round(transT/deltaT);
        T_n         = round(T/deltaT);
        
        % calculate transition probability matrix
        Pij_transT = Pij^transT_n;
        
        % calculate probability to arrive at i at T and derivative
        Pi_T = initProb*Pij^T_n;
        
        % the derivatives are just equal to 0
        % indeed, for Markov_P these are not needed since the delivery time
        % is fixed
        Pij_transT_dot  = zeros(motionModel.indices.nSubPhases,motionModel.indices.nSubPhases);
        Pi_T_dot        = zeros(1,motionModel.indices.nSubPhases);
        
    case 'single'
        
        % allocate the probability matrices and vectors
        Pij_transT  = zeros(motionModel.numPhases,motionModel.numPhases);
        Pi_T        = zeros(1,motionModel.numPhases);
        
        % determine phases at beginning and end of arc
        phase_I = motionModel.lSimulated(abs(motionModel.tSimulated-T) <= 1e-6.*T);
        phase_F = motionModel.lSimulated(abs(motionModel.tSimulated-(T+transT)) <= 1e-6.*(T+transT));
        
        if isempty(phase_I) || isempty(phase_F)
            error('Phases not found!');
        end
        
        % use these phases to construct the probability matrices and vectors
        Pij_transT(phase_I,phase_F) = 1;
        Pi_T(phase_I)               = 1;
        
        % the derivatives are just equal to 0
        Pij_transT_dot  = zeros(motionModel.numPhases,motionModel.numPhases);
        Pi_T_dot        = zeros(1,motionModel.numPhases);
        
    case 'precalculated'
        
        % round the transT and T
        transT  = round(transT,6);
        T       = round(T,6);
        
        % find the appropriate Pij matrix and Pi vector from the
        % precalculated set
        Pij_transT  = motionModel.Pij_transT(:,:,motionModel.transT == transT);
        Pi_T        = motionModel.Pi_arrivalT(:,motionModel.arrivalT == T)';
        
        % the derivatives can just be equal to 0
        % indeed, for precalculated these are not needed since the delivery
        % time is fixed
        Pij_transT_dot  = zeros(motionModel.indices.nSubPhases,motionModel.indices.nSubPhases);
        Pi_T_dot        = zeros(1,motionModel.indices.nSubPhases);
end

end