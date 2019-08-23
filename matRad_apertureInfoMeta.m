function apertureInfo = matRad_apertureInfoMeta(apertureInfo,pln,stf,totalNumOfOptBixels)



% Jacobian matrix to be used in the DAO gradient function
% this tells us the gradient of a particular bixel with respect to an
% element in the apertureVector (aperture weight or leaf position)
% store as a vector for now, convert to sparse matrix later
optBixelFactor = 7;
% For optimized beams: 7 = (1 from weights) + (3 from left leaf positions (I, M, and F)) + (3 from
% right leaf positions (I, M, and F))


if pln.propOpt.runVMAT
    
    tempStruct_beam = apertureInfo.propVMAT.beam;
    tempStruct_jacobT = apertureInfo.propVMAT.jacobT;
    apertureInfo.propVMAT = pln.propOpt.VMAToptions;
    apertureInfo.propVMAT.beam = tempStruct_beam;
    apertureInfo.propVMAT.jacobT = tempStruct_jacobT;
    
    apertureInfo.totalNumOfOptBixels = totalNumOfOptBixels;
    apertureInfo.doseTotalNumOfLeafPairs = sum([apertureInfo.beam(:).numOfActiveLeafPairs]);
    
    % individual size of jacobian vector
    % For interpolated beams: multiply this number times 2 (influenced by the
    % one before and the one after), then add 2 (influenced by the time of the
    % times before and after)
    intBixelFactor = 2*optBixelFactor+2;
    % for the time (probability) gradients
    optBixelFactor = optBixelFactor+apertureInfo.totalNumOfShapes;
    intBixelFactor = intBixelFactor+apertureInfo.totalNumOfShapes;
    
    if apertureInfo.propVMAT.continuousAperture
        apertureInfo.totalNumOfLeafPairs = sum(reshape([apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).doseAngleDAO],2,[]),1)*[apertureInfo.beam([apertureInfo.propVMAT.beam.DAOBeam]).numOfActiveLeafPairs]';
        
        shapeInd = 1;
        for i = 1:numel(apertureInfo.beam)
            if apertureInfo.propVMAT.beam(i).DAOBeam
                
                %apertureInfo.propVMAT.beam(i).initialLeftLeafInd = apertureInfo.beam(i).shape{1}(1).vectorOffset:(apertureInfo.beam(i).shape{1}(1).vectorOffset+dimZ-1);
                %apertureInfo.propVMAT.beam(i).finalLeftLeafInd = (apertureInfo.beam(i).shape{1}(1).vectorOffset+dimZ):(apertureInfo.beam(i).shape{1}(1).vectorOffset+2*dimZ-1);
                
                %apertureInfo.propVMAT.beam(i).initialRightLeafInd = apertureInfo.propVMAT.beam(i).initialLeftLeafInd+apertureInfo.totalNumOfLeafPairs;
                %apertureInfo.propVMAT.beam(i).finalRightLeafInd = apertureInfo.propVMAT.beam(i).finalLeftLeafInd+apertureInfo.totalNumOfLeafPairs;
                
                apertureInfo.propVMAT.beam(i).timeInd = apertureInfo.numPhases*(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)+shapeInd;
                
                %{
                apertureInfo.propVMAT.beam(i).timeInd = repmat(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2+shapeInd,1,apertureInfo.beam(1).numOfActiveLeafPairs);
                
                if apertureInfo.propVMAT.beam(i).timeFac(1) ~= 0
                    % some of the time for this optimized angle "leaks"
                    % into a previous dose angle arc (where we need to
                    % constrain leaf speed)
                    
                    % the initial leaves here act as final leaves for the
                    % previous angle
                    %finalLeftLeafInd = [apertureInfo.propVMAT.beam(i).initialLeftLeafInd apertureInfo.propVMAT.beam(i).finalLeftLeafInd];
                    %finalRightLeafInd = [apertureInfo.propVMAT.beam(i).initialRightLeafInd apertureInfo.propVMAT.beam(i).finalRightLeafInd];
                    
                    % part of the optimized time is used for the leaf speed
                    % calc
                    apertureInfo.propVMAT.beam(i).timeInd = [repmat(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2+shapeInd,1,apertureInfo.beam(1).numOfActiveLeafPairs) apertureInfo.propVMAT.beam(i).timeInd];
                else
                    %finalLeftLeafInd = apertureInfo.propVMAT.beam(i).finalLeftLeafInd;
                    %finalRightLeafInd = apertureInfo.propVMAT.beam(i).finalRightLeafInd;
                end
                
                if apertureInfo.propVMAT.beam(i).timeFac(3) ~= 0
                    % some of the time for this optimized angle "leaks"
                    % into a next dose angle arc (where we need to
                    % constrain leaf speed)
                    
                    %initialLeftLeafInd = [apertureInfo.propVMAT.beam(i).initialLeftLeafInd apertureInfo.propVMAT.beam(i).finalLeftLeafInd];
                    %initialRightLeafInd = [apertureInfo.propVMAT.beam(i).initialRightLeafInd apertureInfo.propVMAT.beam(i).finalRightLeafInd];
                    
                    % part of the optimized time is used for the leaf speed
                    % calc
                    apertureInfo.propVMAT.beam(i).timeInd = [apertureInfo.propVMAT.beam(i).timeInd repmat(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2+shapeInd,1,apertureInfo.beam(1).numOfActiveLeafPairs)];
                else
                    %initialLeftLeafInd = apertureInfo.propVMAT.beam(i).initialLeftLeafInd;
                    %initialRightLeafInd = apertureInfo.propVMAT.beam(i).initialRightLeafInd;
                end
                
                %apertureInfo.propVMAT.beam(i).initialLeftLeafInd = initialLeftLeafInd;
                %apertureInfo.propVMAT.beam(i).initialRightLeafInd = initialRightLeafInd;
                %apertureInfo.propVMAT.beam(i).finalLeftLeafInd = finalLeftLeafInd;
                %apertureInfo.propVMAT.beam(i).finalRightLeafInd = finalRightLeafInd;
                %}
                
                shapeInd = shapeInd+1;
            end
        end
        
    else
        apertureInfo.totalNumOfLeafPairs = totalNumOfLeafPairs;
    end
    
    % fix instances of leaf touching
    apertureInfo = matRad_leafTouching(apertureInfo);
    
    if pln.propOpt.run4D
        apertureInfo = matRad_doDAD(apertureInfo,stf);
        
        % prepare motion model
        apertureInfo.motionModel = matRad_prepModelForOpt(pln.propOpt.prop4D);
        
    else
        
        % use "dummy" motion model, which gives probability 1 for trivial
        % trajectory
        apertureInfo.motionModel.type   = 'Markov';
        apertureInfo.motionModel.qij       = 0;
        apertureInfo.motionModel.initProb  = 1;
        
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
    
    interpGetsTrans = false;
    % prep for gradient offset
    gradOffset      = 0;
    for i = 1:numel(apertureInfo.beam)
        
        maxTime = apertureInfo.propVMAT.beam(i).doseAngleBordersDiff./machine.constraints.gantryRotationSpeed(1);
        [Pij_transT,~,~,~] = matRad_transAndTProb(maxTime,0,apertureInfo.motionModel);
        PIJ_transT = accumarray([apertureInfo.motionModel.indices.subPhase2PosPhase_gridI(:) apertureInfo.motionModel.indices.subPhase2PosPhase_gridJ(:)],Pij_transT(:))./apertureInfo.motionModel.indices.nSubPhasePerPosPhase;
        transMask = PIJ_transT > 0.01;
        
        apertureInfo.propVMAT.beam(i).transMask             = repmat(1:apertureInfo.numPhases,[apertureInfo.numPhases 1]);
        apertureInfo.propVMAT.beam(i).transMask(~transMask) = 0;
        
        if apertureInfo.propVMAT.continuousAperture
            if apertureInfo.propVMAT.beam(i).DAOBeam || interpGetsTrans
                
                apertureInfo.propVMAT.beam(i).leafConstMask             = repmat(1:apertureInfo.numPhases,[apertureInfo.numPhases 1]);
                apertureInfo.propVMAT.beam(i).leafConstMask(~transMask) = 0;
            end
            if apertureInfo.propVMAT.beam(i).DAOBeam
                interpGetsTrans = apertureInfo.propVMAT.beam(i).timeFac(3) ~= 0;
            else
                interpGetsTrans = false;
            end
        end
        
        % individual size of jacobian vector
        if apertureInfo.propVMAT.beam(i).DAOBeam
            apertureInfo.beam(i).bixelJApVec_sz = nnz(~isnan(apertureInfo.beam(i).bixelIndMap)).*optBixelFactor.*apertureInfo.numPhases;
            apertureInfo.beam(i).numUniqueVar   = apertureInfo.beam(i).numOfShapes.*apertureInfo.numPhases.*(1+4.*apertureInfo.beam(i).numOfActiveLeafPairs)+apertureInfo.totalNumOfShapes;%apertureInfo.beam(i).numOfShapes.*(1+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*(apertureInfo.numPhases+1))+apertureInfo.totalNumOfShapes;
            
            % conversion from local opt variable to global opt variable
            apertureInfo.beam(i).local2GlobalVar = zeros(apertureInfo.beam(i).numUniqueVar,1);
            
            for phase = 1:apertureInfo.numPhases
                % weight variables
                apertureInfo.beam(i).local2GlobalVar(phase) = apertureInfo.propVMAT.beam(i).DAOIndex+(phase-1)*apertureInfo.totalNumOfShapes;
                % leaf position variables
                apertureInfo.beam(i).local2GlobalVar(apertureInfo.numPhases+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*(phase-1)+(1:apertureInfo.beam(i).numOfActiveLeafPairs)) ...
                    = apertureInfo.beam(i).shape{phase}(1).vectorOffset(1) + ((1:apertureInfo.beam(i).numOfActiveLeafPairs)-1);
                apertureInfo.beam(i).local2GlobalVar(apertureInfo.numPhases+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*(phase-1)+apertureInfo.beam(i).numOfActiveLeafPairs+(1:apertureInfo.beam(i).numOfActiveLeafPairs)) ...
                    = apertureInfo.beam(i).shape{phase}(1).vectorOffset(2) + ((1:apertureInfo.beam(i).numOfActiveLeafPairs)-1);
                apertureInfo.beam(i).local2GlobalVar(apertureInfo.numPhases+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*(phase-1)+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*apertureInfo.numPhases+(1:apertureInfo.beam(i).numOfActiveLeafPairs)) ...
                    = apertureInfo.beam(i).shape{phase}(1).vectorOffset(1) + apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases + ((1:apertureInfo.beam(i).numOfActiveLeafPairs)-1);
                apertureInfo.beam(i).local2GlobalVar(apertureInfo.numPhases+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*(phase-1)+apertureInfo.beam(i).numOfActiveLeafPairs+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*apertureInfo.numPhases+(1:apertureInfo.beam(i).numOfActiveLeafPairs)) ...
                    = apertureInfo.beam(i).shape{phase}(1).vectorOffset(2) + apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases + ((1:apertureInfo.beam(i).numOfActiveLeafPairs)-1);
                % time variables
                apertureInfo.beam(i).local2GlobalVar(apertureInfo.beam(i).numUniqueVar-apertureInfo.totalNumOfShapes+(1:apertureInfo.totalNumOfShapes)) = ...
                    (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases+(1:apertureInfo.totalNumOfShapes);
            end
            
            if apertureInfo.varOpt
                % keep weight and time variables when calculating gradients for
                % d2
                apertureInfo.beam(i).d2KeepVar = true(apertureInfo.beam(i).numUniqueVar,1);
                %apertureInfo.beam(i).d2KeepVar((apertureInfo.beam(i).numOfShapes.*apertureInfo.numPhases+1):end) = false;
                %apertureInfo.beam(i).d2KeepVar((apertureInfo.beam(i).numOfShapes.*apertureInfo.numPhases+1):(end-apertureInfo.totalNumOfShapes)) = false;
            end
            
        else
            apertureInfo.beam(i).bixelJApVec_sz = nnz(~isnan(apertureInfo.beam(i).bixelIndMap)).*intBixelFactor.*apertureInfo.numPhases;
            apertureInfo.beam(i).numUniqueVar   = apertureInfo.numPhases.*(2+4.*apertureInfo.beam(i).numOfActiveLeafPairs)+apertureInfo.totalNumOfShapes;%(2+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*(apertureInfo.numPhases+1))+apertureInfo.totalNumOfShapes;
            
            % conversion from local opt variable to global opt variable
            apertureInfo.beam(i).local2GlobalVar = zeros(apertureInfo.beam(i).numUniqueVar,1);
            
            for phase = 1:apertureInfo.numPhases
                % weight variables
                apertureInfo.beam(i).local2GlobalVar(2.*(phase-1)+(1:2)) = [apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).DAOIndex apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).DAOIndex]+(phase-1)*apertureInfo.totalNumOfShapes;
                % leaf position variables
                apertureInfo.beam(i).local2GlobalVar(2.*apertureInfo.numPhases+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*(phase-1)+(1:apertureInfo.beam(i).numOfActiveLeafPairs)) ...
                    = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(1).vectorOffset(2) + ((1:apertureInfo.beam(i).numOfActiveLeafPairs)-1);
                apertureInfo.beam(i).local2GlobalVar(2.*apertureInfo.numPhases+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*(phase-1)+apertureInfo.beam(i).numOfActiveLeafPairs+(1:apertureInfo.beam(i).numOfActiveLeafPairs)) ...
                    = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(1).vectorOffset(1) + ((1:apertureInfo.beam(i).numOfActiveLeafPairs)-1);
                apertureInfo.beam(i).local2GlobalVar(2.*apertureInfo.numPhases+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*(phase-1)+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*apertureInfo.numPhases+(1:apertureInfo.beam(i).numOfActiveLeafPairs)) ...
                    = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(1).vectorOffset(2) + apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases + ((1:apertureInfo.beam(i).numOfActiveLeafPairs)-1);
                apertureInfo.beam(i).local2GlobalVar(2.*apertureInfo.numPhases+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*(phase-1)+apertureInfo.beam(i).numOfActiveLeafPairs+2.*apertureInfo.beam(i).numOfActiveLeafPairs.*apertureInfo.numPhases+(1:apertureInfo.beam(i).numOfActiveLeafPairs)) ...
                    = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(1).vectorOffset(1) + apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases + ((1:apertureInfo.beam(i).numOfActiveLeafPairs)-1);
                % time variables
                apertureInfo.beam(i).local2GlobalVar(apertureInfo.beam(i).numUniqueVar-apertureInfo.totalNumOfShapes+(1:apertureInfo.totalNumOfShapes)) = ...
                    (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases+(1:apertureInfo.totalNumOfShapes);
            end
            
            if apertureInfo.varOpt
                % keep weight and time variables when calculating gradients for
                % d2
                apertureInfo.beam(i).d2KeepVar = true(apertureInfo.beam(i).numUniqueVar,1);
                %apertureInfo.beam(i).d2KeepVar((2.*apertureInfo.numPhases+1):end) = false;
                %apertureInfo.beam(i).d2KeepVar((2.*apertureInfo.numPhases+1):(end-apertureInfo.totalNumOfShapes)) = false;
            end
        end
        
        if apertureInfo.varOpt
            
            apertureInfo.beam(i).d2KeepVar = find(apertureInfo.beam(i).d2KeepVar);
            apertureInfo.beam(i).numKeepVar = numel(apertureInfo.beam(i).d2KeepVar);
            
            apertureInfo.beam(i).gradOffset = gradOffset;
            gradOffset = gradOffset+apertureInfo.beam(i).numKeepVar.*apertureInfo.numPhases.^2;
        end
    end
    
    if apertureInfo.propVMAT.continuousAperture
        
        % count number of transitions
        apertureInfo.propVMAT.numLeafSpeedConstraint      = nnz([apertureInfo.propVMAT.beam.leafConstMask]);
        apertureInfo.propVMAT.numLeafSpeedConstraintDAO   = nnz([apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).leafConstMask]);
    end
    
else
    apertureInfo.totalNumOfLeafPairs = sum([apertureInfo.beam.numOfShapes]*[apertureInfo.beam.numOfActiveLeafPairs]');
    apertureInfo.doseTotalNumOfLeafPairs = apertureInfo.totalNumOfLeafPairs;
    
    % individual size of jacobian vector
    % For interpolated beams: multiply this number times 2 (influenced by the
    % one before and the one after)
    intBixelFactor = 2*optBixelFactor;
    
    for i = 1:numel(apertureInfo.beam)
        if apertureInfo.propVMAT.beam(i).DAOBeam
            apertureInfo.beam(i).bixelJApVec_sz = nnz(~isnan(apertureInfo.beam(i).bixelIndMap)).*optBixelFactor.*apertureInfo.beam(i).numOfShapes.*apertureInfo.numPhases;
        else
            apertureInfo.beam(i).bixelJApVec_sz = nnz(~isnan(apertureInfo.beam(i).bixelIndMap)).*intBixelFactor.*apertureInfo.beam(i).numOfShapes.*apertureInfo.numPhases;
        end
    end
end

% global size of jacobian vector
apertureInfo.bixelJApVec_sz = sum([apertureInfo.beam.bixelJApVec_sz]);
% apertureInfo.bixelJApVec_sz = (apertureInfo.totalNumOfOptBixels*optBixelFactor+(apertureInfo.totalNumOfBixels-apertureInfo.totalNumOfOptBixels)*intBixelFactor)*apertureInfo.numPhases;

end

