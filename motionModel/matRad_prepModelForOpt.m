function motionModel = matRad_prepModelForOpt(pln,stf,apertureInfo)

% check if motion model already exists
if isfield(pln.propOpt.prop4D,'motionModel')
    
    % if we're doing a fixed gantry speed, we can pre-calculate all of the
    % probabilities to save time
    if pln.propOpt.VMAToptions.fixedGantrySpeed
        
        % fixed gantry speed, use transition probability matrix P to
        % precalculate probabilities
        
        % use model if it exists (and we're not doing a fixed gantry speed)
        motionModel        = pln.propOpt.prop4D.motionModel;
        motionModel.type   = 'Markov_P';
        
        % strip model of the "out of bounds" phases
        %motionModel = matRad_stripMarkovOOB(motionModel);
        
        % get gridded maps
        [motionModel.indices.subPhase2PosPhase_gridJ, motionModel.indices.subPhase2PosPhase_gridI] = meshgrid(motionModel.indices.subPhase2PosPhase);
        
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
        
        
        % calculate the fixed gantry speed
        % this should be equal to the speed calculated in matRad_arcSequencing
        %gantryRot = (pln.propOpt.VMAToptions.finishingAngle-pln.propOpt.VMAToptions.startingAngle)./pln.propOpt.VMAToptions.deliveryTime;
        gantryRot = apertureInfo.beam(1).gantryRot;
        
        % loop through the stf structure to find all possible transition and
        % arrival times
        arrivalT        = zeros(numel(stf),1);
        transT_i1       = zeros(numel(stf),1);
        transT_i1DAO    = zeros(numel(pln.propStf.DAOGantryAngles),1);
        transT_i1i2     = zeros(numel(stf).*(numel(stf)-1)./2,1);
        i1DAO           = 1;
        i1i2            = 1;
        
        % loop through beams i1
        for i1 = 1:numel(stf)
            
            % calculate arrival time for i1
            arrivalT(i1) = (stf(i1).propVMAT.fluAngleBorders(1)-pln.propOpt.VMAToptions.startingAngle)./gantryRot;
            
            % calculate transitiion time for i1
            transT_i1(i1) = stf(i1).propVMAT.fluAngleBordersDiff./gantryRot;
            
            if stf(i1).propVMAT.DAOBeam
                transT_i1DAO(i1DAO) = stf(i1).propVMAT.DAOAngleBordersDiff./gantryRot;
                
                i1DAO = i1DAO+1;
            end
            
            % loop through beams i2
            for i2 = (i1+1):numel(stf)
                
                % calculate transition time from i1 to i2
                transT_i1i2(i1i2) = (stf(i2).propVMAT.fluAngleBorders(1)-stf(i1).propVMAT.fluAngleBorders(2))./gantryRot;
                
                i1i2 = i1i2+1;
            end
        end
        
        % get unique (rounded) values for both types of transition times
        % also round arrival times
        transT_unique   = unique(round([transT_i1; transT_i1DAO; transT_i1i2],6));
        arrivalT        = round(arrivalT,6);
        
        % put new motionModel stuff in:
        % 1. arrival times and corresponding probability vectors
        % 2. transition times and corresponding probability matrices
        motionModel.arrivalT    = arrivalT;
        motionModel.Pi_arrivalT = zeros(motionModel.indices.nSubPhases,numel(arrivalT));
        motionModel.transT      = transT_unique;
        motionModel.Pij_transT  = zeros(motionModel.indices.nSubPhases,motionModel.indices.nSubPhases,numel(transT_unique));
        
        % loop through all arrival time values to form the probability
        % vectors
        for i = 1:numel(arrivalT)
            
            % calculate probability vector for current arrival time
            [~,~,Pi_arrivalT,~] = matRad_transAndTProb(0,arrivalT(i),motionModel);
            
            % insert into motionModel structure
            motionModel.Pi_arrivalT(:,i) = Pi_arrivalT';
        end
        
        % loop through all transition times to form the probability
        % matrices
        for i = 1:numel(transT_unique)
            
            % calculate probability vector for current arrival time
            [Pij_transT,~,~,~] = matRad_transAndTProb(transT_unique(i),0,motionModel);
            
            % insert into motionModel structure
            motionModel.Pij_transT(:,:,i) = Pij_transT;
        end
        
        % change the motionModel type to precalculated
        motionModel.type = 'precalculated';
        
    else
        % variable gantry speed, use transition rate matrix q
        % not currently supported
        
        % use model if it exists (and we're not doing a fixed gantry speed)
        motionModel        = pln.propOpt.prop4D.motionModel;
        motionModel.type   = 'Markov_Q';
        
        % strip model of the "out of bounds" phases
        %motionModel = matRad_stripMarkovOOB(motionModel);
        
        % get gridded maps
        [motionModel.indices.subPhase2PosPhase_gridJ, motionModel.indices.subPhase2PosPhase_gridI] = meshgrid(motionModel.indices.subPhase2PosPhase);
        
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
    end
    
else
    
    % if it doesn't, throw an error
    error('No motion model found in structure pln.propOpt.prop4d');
end

end

