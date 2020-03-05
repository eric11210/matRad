function apertureInfo = matRad_apertureInfoMeta(apertureInfo,pln,stf,doLeafStuff)

if nargin < 4
    doLeafStuff = true;
end


% Jacobian matrix to be used in the DAO gradient function
% this tells us the gradient of a particular bixel with respect to an
% element in the apertureVector (aperture weight or leaf position)
% store as a vector for now, convert to sparse matrix later

% for individual beam size of jacobian vector
optNumVarMult = 9;
% For optimized beams: 9 = (1 from weights) + (4 from left leaf positions)
% + (4 from right leaf positions)
if pln.propOpt.VMAToptions.fixedGantrySpeed
    intNumVarMult = 10;
else
    intNumVarMult = 12;
end
% For interpolated beams: 12 = (2 from weights last/next) + (4 from left
% leaf positions) + (4 from right leaf positions) + (2 from times 
% last/next)

% initialize number of fixels
lastFixelIndOffset = 0;
nextFixelIndOffset = 0;

if pln.propOpt.runVMAT
    
    tempStruct_beam = apertureInfo.propVMAT.beam;
    tempStruct_jacobT = apertureInfo.propVMAT.jacobT;
    apertureInfo.propVMAT = pln.propOpt.VMAToptions;
    apertureInfo.propVMAT.beam = tempStruct_beam;
    apertureInfo.propVMAT.jacobT = tempStruct_jacobT;
    
    apertureInfo.doseTotalNumOfLeafPairs = sum([apertureInfo.beam(:).numOfActiveLeafPairs]);
    
    if apertureInfo.propVMAT.continuousAperture
        apertureInfo.totalNumOfLeafPairs = sum(reshape([apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).doseAngleDAO],2,[]),1)*[apertureInfo.beam([apertureInfo.propVMAT.beam.DAOBeam]).numOfActiveLeafPairs]';
        
        shapeInd = 1;
        for i = 1:numel(apertureInfo.beam)
            
            % determine last/next bixelIndMaps
            apertureInfo.beam(i).lastBixelIndMap = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDoseIndex).bixelIndMap;
            apertureInfo.beam(i).nextBixelIndMap = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDoseIndex).bixelIndMap;
            
            % determine last/next number of bixels
            apertureInfo.beam(i).lastNumBixels = nnz(~isnan(apertureInfo.beam(i).lastBixelIndMap));
            apertureInfo.beam(i).nextNumBixels = nnz(~isnan(apertureInfo.beam(i).nextBixelIndMap));
            
            % determine last/next fixelIndMaps
            apertureInfo.beam(i).lastFixelIndMap                                                = NaN(size(apertureInfo.beam(i).lastBixelIndMap));
            apertureInfo.beam(i).lastFixelIndMap(~isnan(apertureInfo.beam(i).lastBixelIndMap))  = (1:apertureInfo.beam(i).lastNumBixels)+lastFixelIndOffset;
            apertureInfo.beam(i).nextFixelIndMap                                                = NaN(size(apertureInfo.beam(i).nextBixelIndMap));
            apertureInfo.beam(i).nextFixelIndMap(~isnan(apertureInfo.beam(i).nextBixelIndMap))  = (1:apertureInfo.beam(i).nextNumBixels)+nextFixelIndOffset;
            
            % update offsets
            lastFixelIndOffset = lastFixelIndOffset+apertureInfo.beam(i).lastNumBixels;
            nextFixelIndOffset = nextFixelIndOffset+apertureInfo.beam(i).nextNumBixels;
            
            % determine effective number of bixels
            apertureInfo.beam(i).effNumBixels_lastDose_arcI = any(apertureInfo.propVMAT.beam(i).fracToLastDose_arcI).*apertureInfo.beam(i).lastNumBixels;
            apertureInfo.beam(i).effNumBixels_lastDose_arcF = any(apertureInfo.propVMAT.beam(i).fracToLastDose_arcF).*apertureInfo.beam(i).lastNumBixels;
            apertureInfo.beam(i).effNumBixels_nextDose_arcI = any(apertureInfo.propVMAT.beam(i).fracToNextDose_arcI).*apertureInfo.beam(i).nextNumBixels;
            apertureInfo.beam(i).effNumBixels_nextDose_arcF = any(apertureInfo.propVMAT.beam(i).fracToNextDose_arcF).*apertureInfo.beam(i).nextNumBixels;
            apertureInfo.beam(i).effNumBixels_lastDose      = max(apertureInfo.beam(i).effNumBixels_lastDose_arcI,apertureInfo.beam(i).effNumBixels_lastDose_arcF);
            apertureInfo.beam(i).effNumBixels_nextDose      = max(apertureInfo.beam(i).effNumBixels_nextDose_arcI,apertureInfo.beam(i).effNumBixels_nextDose_arcF);
            
            if apertureInfo.propVMAT.beam(i).DAOBeam
                
                apertureInfo.propVMAT.beam(i).timeInd = apertureInfo.numPhases*(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)+shapeInd;
                
                shapeInd = shapeInd+1;
            end
        end
        
        % record number of fixels
        apertureInfo.totalNumOfLastFixels = lastFixelIndOffset;
        apertureInfo.totalNumOfNextFixels = nextFixelIndOffset;
        
    else
        apertureInfo.totalNumOfLeafPairs = totalNumOfLeafPairs;
    end
    
    if doLeafStuff
        % fix instances of leaf touching
        apertureInfo = matRad_leafTouching(apertureInfo);
    end
    
    if pln.propOpt.run4D
        
        if doLeafStuff
            apertureInfo = matRad_doDAD(apertureInfo,stf);
        end
        
        % prepare motion model
        apertureInfo.motionModel = matRad_prepModelForOpt(pln,stf,apertureInfo);
        
    else
        
        % use "dummy" motion model, which gives probability 1 for trivial
        % trajectory
        apertureInfo.motionModel.type       = 'Markov_Q';
        apertureInfo.motionModel.qij        = 0;
        apertureInfo.motionModel.initProb   = 1;
        
        % diagonalize matrix
        [apertureInfo.motionModel.qij_V,apertureInfo.motionModel.qij_D] = eig(apertureInfo.motionModel.qij);
        
        apertureInfo.motionModel.indices.subPhase2PosPhase_gridI    = 1;
        apertureInfo.motionModel.indices.subPhase2PosPhase_gridJ    = 1;
        apertureInfo.motionModel.indices.nSubPhasePerPosPhase       = 1;
        apertureInfo.motionModel.indices.nSubPhases                 = 1;
    end
    
    fileName = apertureInfo.propVMAT.machineConstraintFile;
    try
        load(fileName,'machine');
    catch
        error(['Could not find the following machine file: ' fileName ]);
    end
    
    % set the minimum gantry rotation speed
    if apertureInfo.propVMAT.fixedGantrySpeed
        % this should correspond to the speed calculated in
        % matRad_arcSequencing
        %minGantryRot = (pln.propOpt.VMAToptions.finishingAngle-pln.propOpt.VMAToptions.startingAngle)./pln.propOpt.VMAToptions.deliveryTime;
        minGantryRot = apertureInfo.beam(1).gantryRot;
    else
        % just use the machine speed
        minGantryRot = machine.constraints.gantryRotationSpeed(1);
    end
    
    % prep for gradient offset
    gradOffset      = 0;
    for i = 1:numel(apertureInfo.beam)
        
        % determine possible phase transitions between fluence angles
        maxTime_flu = apertureInfo.propVMAT.beam(i).fluAngleBordersDiff./minGantryRot;
        [Pij_transT_flu,~,~,~] = matRad_transAndTProb(maxTime_flu,0,apertureInfo.motionModel);
        PIJ_transT_flu = accumarray([apertureInfo.motionModel.indices.subPhase2PosPhase_gridI(:) apertureInfo.motionModel.indices.subPhase2PosPhase_gridJ(:)],Pij_transT_flu(:))./apertureInfo.motionModel.indices.nSubPhasePerPosPhase;
        transMask_flu = PIJ_transT_flu > 0;
        
        apertureInfo.propVMAT.beam(i).transMask                 = repmat(1:apertureInfo.numPhases,[apertureInfo.numPhases 1]);
        apertureInfo.propVMAT.beam(i).transMask(~transMask_flu) = 0;
        
        % calculate individual size of jacobian vector
        % construct conversion from local to global optimization variables
        if apertureInfo.propVMAT.beam(i).DAOBeam
            apertureInfo.beam(i).bixelJApVecLastDose_sz = apertureInfo.beam(i).effNumBixels_lastDose.*(optNumVarMult+apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).DAOIndex);
            apertureInfo.beam(i).bixelJApVecNextDose_sz = apertureInfo.beam(i).effNumBixels_nextDose.*(optNumVarMult+apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).DAOIndex);
            if pln.propOpt.VMAToptions.fixedGantrySpeed
                apertureInfo.beam(i).numFullVar             = apertureInfo.beam(i).numOfShapes.*apertureInfo.numPhases.*(1+8.*apertureInfo.beam(i).numOfActiveLeafPairs);%apertureInfo.beam(i).numOfShapes.*(1+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*(apertureInfo.numPhases+1))+apertureInfo.totalNumOfShapes;
            else
                apertureInfo.beam(i).numFullVar             = apertureInfo.beam(i).numOfShapes.*apertureInfo.numPhases.*(1+8.*apertureInfo.beam(i).numOfActiveLeafPairs)+apertureInfo.totalNumOfShapes;%apertureInfo.beam(i).numOfShapes.*(1+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*(apertureInfo.numPhases+1))+apertureInfo.totalNumOfShapes;
            end
            % conversion from local opt variable to global opt variable
            apertureInfo.beam(i).local2GlobalVar = zeros(apertureInfo.beam(i).numFullVar,1);
            
            for phase = 1:apertureInfo.numPhases
                % weight variables
                apertureInfo.beam(i).local2GlobalVar(phase) = apertureInfo.propVMAT.beam(i).DAOIndex+(phase-1)*apertureInfo.totalNumOfShapes;
            end
            
            % calculate weight offset for local2GlobalVar
            weightOffset = apertureInfo.numPhases;
        else
            apertureInfo.beam(i).bixelJApVecLastDose_sz = apertureInfo.beam(i).effNumBixels_lastDose.*(intNumVarMult+apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).DAOIndex);
            apertureInfo.beam(i).bixelJApVecNextDose_sz = apertureInfo.beam(i).effNumBixels_nextDose.*(intNumVarMult+apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).DAOIndex);
            if pln.propOpt.VMAToptions.fixedGantrySpeed
                apertureInfo.beam(i).numFullVar             = apertureInfo.numPhases.*(2+8.*apertureInfo.beam(i).numOfActiveLeafPairs);%(2+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*(apertureInfo.numPhases+1))+apertureInfo.totalNumOfShapes;
            else
                apertureInfo.beam(i).numFullVar             = apertureInfo.numPhases.*(2+8.*apertureInfo.beam(i).numOfActiveLeafPairs)+apertureInfo.totalNumOfShapes;%(2+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*(apertureInfo.numPhases+1))+apertureInfo.totalNumOfShapes;
            end
            
            % conversion from local opt variable to global opt variable
            apertureInfo.beam(i).local2GlobalVar = zeros(apertureInfo.beam(i).numFullVar,1);
            
            for phase = 1:apertureInfo.numPhases
                % weight variables
                apertureInfo.beam(i).local2GlobalVar(2.*(phase-1)+(1:2)) = [apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).DAOIndex apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).DAOIndex]+(phase-1)*apertureInfo.totalNumOfShapes;
            end
            
            % calculate weight offset for local2GlobalVar
            weightOffset = 2.*apertureInfo.numPhases;
        end
        
        % leaf position variables
        for phase = 1:apertureInfo.numPhases
            % L_lastDAOI
            apertureInfo.beam(i).local2GlobalVar(weightOffset+4.*apertureInfo.beam(i).numOfActiveLeafPairs.*(phase-1)+0.*apertureInfo.beam(i).numOfActiveLeafPairs+(1:apertureInfo.beam(i).numOfActiveLeafPairs)) ...
                = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(1).vectorOffset(1) + ((1:apertureInfo.beam(i).numOfActiveLeafPairs)-1);
            % L_lastDAOF
            apertureInfo.beam(i).local2GlobalVar(weightOffset+4.*apertureInfo.beam(i).numOfActiveLeafPairs.*(phase-1)+1.*apertureInfo.beam(i).numOfActiveLeafPairs+(1:apertureInfo.beam(i).numOfActiveLeafPairs)) ...
                = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(1).vectorOffset(2) + ((1:apertureInfo.beam(i).numOfActiveLeafPairs)-1);
            % L_nextDAOI
            apertureInfo.beam(i).local2GlobalVar(weightOffset+4.*apertureInfo.beam(i).numOfActiveLeafPairs.*(phase-1)+2.*apertureInfo.beam(i).numOfActiveLeafPairs+(1:apertureInfo.beam(i).numOfActiveLeafPairs)) ...
                = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(1).vectorOffset(1) + ((1:apertureInfo.beam(i).numOfActiveLeafPairs)-1);
            % L_nextDAOF
            apertureInfo.beam(i).local2GlobalVar(weightOffset+4.*apertureInfo.beam(i).numOfActiveLeafPairs.*(phase-1)+3.*apertureInfo.beam(i).numOfActiveLeafPairs+(1:apertureInfo.beam(i).numOfActiveLeafPairs)) ...
                = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(1).vectorOffset(2) + ((1:apertureInfo.beam(i).numOfActiveLeafPairs)-1);
            
            % R_lastDAOI
            apertureInfo.beam(i).local2GlobalVar(weightOffset+4.*apertureInfo.beam(i).numOfActiveLeafPairs.*(phase-1)+0.*apertureInfo.beam(i).numOfActiveLeafPairs+4.*apertureInfo.beam(i).numOfActiveLeafPairs.*apertureInfo.numPhases+(1:apertureInfo.beam(i).numOfActiveLeafPairs)) ...
                = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(1).vectorOffset(1) + apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases + ((1:apertureInfo.beam(i).numOfActiveLeafPairs)-1);
            % R_lastDAOF
            apertureInfo.beam(i).local2GlobalVar(weightOffset+4.*apertureInfo.beam(i).numOfActiveLeafPairs.*(phase-1)+1.*apertureInfo.beam(i).numOfActiveLeafPairs+4.*apertureInfo.beam(i).numOfActiveLeafPairs.*apertureInfo.numPhases+(1:apertureInfo.beam(i).numOfActiveLeafPairs)) ...
                = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(1).vectorOffset(2) + apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases + ((1:apertureInfo.beam(i).numOfActiveLeafPairs)-1);
            % R_nextDAOI
            apertureInfo.beam(i).local2GlobalVar(weightOffset+4.*apertureInfo.beam(i).numOfActiveLeafPairs.*(phase-1)+2.*apertureInfo.beam(i).numOfActiveLeafPairs+4.*apertureInfo.beam(i).numOfActiveLeafPairs.*apertureInfo.numPhases+(1:apertureInfo.beam(i).numOfActiveLeafPairs)) ...
                = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(1).vectorOffset(1) + apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases + ((1:apertureInfo.beam(i).numOfActiveLeafPairs)-1);
            % R_nextDAOF
            apertureInfo.beam(i).local2GlobalVar(weightOffset+4.*apertureInfo.beam(i).numOfActiveLeafPairs.*(phase-1)+3.*apertureInfo.beam(i).numOfActiveLeafPairs+4.*apertureInfo.beam(i).numOfActiveLeafPairs.*apertureInfo.numPhases+(1:apertureInfo.beam(i).numOfActiveLeafPairs)) ...
                = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(1).vectorOffset(2) + apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases + ((1:apertureInfo.beam(i).numOfActiveLeafPairs)-1);
        end
        
        if ~pln.propOpt.VMAToptions.fixedGantrySpeed
            % time variables
            apertureInfo.beam(i).local2GlobalVar(apertureInfo.beam(i).numFullVar-apertureInfo.totalNumOfShapes+(1:apertureInfo.totalNumOfShapes)) = ...
                (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases+(1:apertureInfo.totalNumOfShapes);
        end
        
        % only keep unique local2GlobalVar variables
        [apertureInfo.beam(i).local2GlobalVar,~,apertureInfo.beam(i).full2UniqueLocalVar]   = unique(apertureInfo.beam(i).local2GlobalVar);
        apertureInfo.beam(i).numUniqueVar                                                   = numel(apertureInfo.beam(i).local2GlobalVar);
        
        if apertureInfo.varOpt
            % decide which variables to keep when calculating gradients for
            % d2
            % keep all
            apertureInfo.beam(i).d2KeepVar = true(apertureInfo.beam(i).numUniqueVar,1);
            % keep weight and time variables when calculating gradients for
            % d2
            %apertureInfo.beam(i).d2KeepVar((apertureInfo.beam(i).numOfShapes.*apertureInfo.numPhases+1):end) = false;
            %apertureInfo.beam(i).d2KeepVar((apertureInfo.beam(i).numOfShapes.*apertureInfo.numPhases+1):(end-apertureInfo.totalNumOfShapes)) = false;
            
            apertureInfo.beam(i).d2KeepVar = find(apertureInfo.beam(i).d2KeepVar);
            apertureInfo.beam(i).numKeepVar = numel(apertureInfo.beam(i).d2KeepVar);
            
            apertureInfo.beam(i).gradOffset = gradOffset;
            gradOffset = gradOffset+apertureInfo.beam(i).numKeepVar.*apertureInfo.numPhases.^2;
        end
    end
    
    if apertureInfo.propVMAT.continuousAperture
        
        % count number of transitions
        apertureInfo.propVMAT.numLeafSpeedConstraint      = nnz([apertureInfo.propVMAT.beam.transMask]);
        apertureInfo.propVMAT.numLeafSpeedConstraintDAO   = nnz([apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).transMask]);
    end
    
else
    apertureInfo.totalNumOfLeafPairs = sum([apertureInfo.beam.numOfShapes]*[apertureInfo.beam.numOfActiveLeafPairs]');
    apertureInfo.doseTotalNumOfLeafPairs = apertureInfo.totalNumOfLeafPairs;
    
    % individual size of jacobian vector
    % For interpolated beams: multiply this number times 2 (influenced by the
    % one before and the one after)
    intBixelFactor = 2*optNumVarMult;
    
    for i = 1:numel(apertureInfo.beam)
        if apertureInfo.propVMAT.beam(i).DAOBeam
            apertureInfo.beam(i).bixelJApVec_sz = nnz(~isnan(apertureInfo.beam(i).bixelIndMap)).*optNumVarMult.*apertureInfo.beam(i).numOfShapes.*apertureInfo.numPhases;
        else
            apertureInfo.beam(i).bixelJApVec_sz = nnz(~isnan(apertureInfo.beam(i).bixelIndMap)).*intBixelFactor.*apertureInfo.beam(i).numOfShapes.*apertureInfo.numPhases;
        end
    end
end

end

