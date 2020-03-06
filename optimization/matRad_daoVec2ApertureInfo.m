function updatedInfo = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfoVect)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to translate the vector representation of the aperture
% shape and weight into an aperture info struct. At the same time, the
% updated bixel weight vector w is computed and a vector listing the
% correspondence between leaf tips and bixel indices for gradient
% calculation
%
% call
%   [updatedInfo,w,indVect] = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfoVect)
%
% input
%   apertureInfo:     aperture shape info struct
%   apertureInfoVect: aperture weights and shapes parameterized as vector
%
% output
%   updatedInfo: updated aperture shape info struct according to apertureInfoVect
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function to update the apertureInfo struct after the each iteraton of the
% optimization

% initializing variables
updatedInfo = apertureInfo;

updatedInfo.apertureVector = apertureInfoVect;


%% preliminary loop

% we'll always need to do a loop over all beams and phases to clear the
% shapeMap, and sumGradSq _weight and _leaf

% if we're doing VMAT, also take this time to update MU, MU rates, gantry
% rotation speeds, etc.

% loop over all beams
for i = 1:numel(updatedInfo.beam)
    
    % loop over all shapes
    if updatedInfo.runVMAT
        numOfShapes = 1;
    else
        numOfShapes = updatedInfo.beam(i).numOfShapes;
    end
    
    % loop over all phases
    for phase = 1:updatedInfo.numPhases
        
        for j = 1:numOfShapes
            
            % clear shapeMap, and sumGradSq _weight and _leaf
            updatedInfo.beam(i).shape{phase}(j).shapeMap            = zeros(size(updatedInfo.beam(i).lastBixelIndMap));
            updatedInfo.beam(i).shape{phase}(j).sumGradSq_weight    = 0;
            updatedInfo.beam(i).shape{phase}(j).sumGradSq_leaf      = 0;
            
            if updatedInfo.runVMAT && updatedInfo.propVMAT.beam(i).DAOBeam
                
                % update the shape weight
                % rescale the weight from the vector using the previous
                % iteration scaling factor
                updatedInfo.beam(i).shape{phase}(j).weight = apertureInfoVect(updatedInfo.beam(i).shape{phase}(j).weightOffset)./updatedInfo.beam(i).shape{phase}(j).jacobiScale;
                
                updatedInfo.beam(i).shape{phase}(j).MU = updatedInfo.beam(i).shape{phase}(j).weight*updatedInfo.weightToMU;
                if phase == 1
                    if ~updatedInfo.propVMAT.fixedGantrySpeed
                        updatedInfo.beam(i).time        = apertureInfoVect((updatedInfo.totalNumOfShapes+updatedInfo.totalNumOfLeafPairs*2)*updatedInfo.numPhases+updatedInfo.propVMAT.beam(i).DAOIndex)*updatedInfo.propVMAT.beam(i).timeFacCurr;
                        updatedInfo.beam(i).gantryRot   = updatedInfo.propVMAT.beam(i).fluAngleBordersDiff/updatedInfo.beam(i).time;
                    end
                end
                updatedInfo.beam(i).shape{phase}(j).MURate = updatedInfo.beam(i).shape{phase}(j).MU./updatedInfo.beam(i).time;
            end
        end
    end
end

%% bixel/fixel weight and jacobian prep

iBitStr = sprintf('uint%d',2^ceil(log2(log2(numel(apertureInfoVect)))));
jBitStr = sprintf('uint%d',2^ceil(log2(log2(updatedInfo.totalNumOfLastFixels))));

% weights
lastFixWeightBase       = cell(updatedInfo.numPhases.^2,1);
lastFixWeightBase(:)    = {zeros(updatedInfo.totalNumOfLastFixels,1)};
nextFixWeightBase       = cell(updatedInfo.numPhases.^2,1);
nextFixWeightBase(:)    = {zeros(updatedInfo.totalNumOfNextFixels,1)};
% make last/next arcI and arcF structs
updatedInfo.arcI.lastDose.fixelWeights = lastFixWeightBase;
updatedInfo.arcI.nextDose.fixelWeights = nextFixWeightBase;
updatedInfo.arcF.lastDose.fixelWeights = lastFixWeightBase;
updatedInfo.arcF.nextDose.fixelWeights = nextFixWeightBase;

% make real bixel weight struct
updatedInfo.bixelWeights        = cell(updatedInfo.numPhases,1);
updatedInfo.bixelWeights(:)     = {zeros(updatedInfo.totalNumOfBixels,1)};

% weights
bixWeightAndGradBase_lastDose_angle.w = cell(updatedInfo.numPhases.^2,1);
bixWeightAndGradBase_nextDose_angle.w = cell(updatedInfo.numPhases.^2,1);
% jacobian
bixWeightAndGradBase_lastDose_angle.bixelJApVec_vec = cell(updatedInfo.numPhases.^2,1);
bixWeightAndGradBase_nextDose_angle.bixelJApVec_vec = cell(updatedInfo.numPhases.^2,1);
% vector indices
bixWeightAndGradBase_lastDose_angle.bixelJApVec_i = cell(updatedInfo.numPhases.^2,1);
bixWeightAndGradBase_nextDose_angle.bixelJApVec_i = cell(updatedInfo.numPhases.^2,1);
% bixel indices
bixWeightAndGradBase_lastDose_angle.bixelJApVec_j = cell(updatedInfo.numPhases.^2,1);
bixWeightAndGradBase_nextDose_angle.bixelJApVec_j = cell(updatedInfo.numPhases.^2,1);
% offset
bixWeightAndGradBase_lastDose_angle.bixelJApVec_offset = cell(updatedInfo.numPhases.^2,1);
bixWeightAndGradBase_nextDose_angle.bixelJApVec_offset = cell(updatedInfo.numPhases.^2,1);
% make last/next arcI and arcF structs
updatedInfo.arcI.lastDose.bixelJApVec = cell(numel(updatedInfo.beam).*updatedInfo.numPhases.^2,1);
updatedInfo.arcI.nextDose.bixelJApVec = cell(numel(updatedInfo.beam).*updatedInfo.numPhases.^2,1);
updatedInfo.arcF.lastDose.bixelJApVec = cell(numel(updatedInfo.beam).*updatedInfo.numPhases.^2,1);
updatedInfo.arcF.nextDose.bixelJApVec = cell(numel(updatedInfo.beam).*updatedInfo.numPhases.^2,1);

% make real bixel jacobian struct
updatedInfo.bixelJApVec = cell(numel(updatedInfo.beam).*updatedInfo.numPhases,1);

% probabilities
updatedInfo.probI_IJ        = cell(numel(updatedInfo.beam),1);
updatedInfo.probI_Ij        = cell(numel(updatedInfo.beam),1);
updatedInfo.probF_KL        = cell(numel(updatedInfo.beam),1);
updatedInfo.probF_kL        = cell(numel(updatedInfo.beam),1);
updatedInfo.probIGrad_IJ    = cell(numel(updatedInfo.beam),1);
updatedInfo.probIGrad_Ij    = cell(numel(updatedInfo.beam),1);
updatedInfo.probFGrad_KL    = cell(numel(updatedInfo.beam),1);
updatedInfo.probFGrad_kL    = cell(numel(updatedInfo.beam),1);

%% update the shapeMaps
% here the new colimator positions are used to create new shapeMaps that
% now include decimal values instead of binary

% loop over all beams
for i = 1:numel(updatedInfo.beam)
    
    % get dimensions of 2d matrices that store shape/bixel information
    n = updatedInfo.beam(i).numOfActiveLeafPairs;
    
    % loop over all shapes
    if updatedInfo.runVMAT
        numOfShapes = 1;
    else
        numOfShapes = updatedInfo.beam(i).numOfShapes;
    end
    
    % loop over all (initial) phases
    for phase = 1:updatedInfo.numPhases
        
        for j = 1:numOfShapes
            
            if ~updatedInfo.runVMAT || updatedInfo.propVMAT.beam(i).DAOBeam
                % either this is not VMAT, or if it is VMAT, this is a DAO beam
                
                % update the shape weight
                updatedInfo.beam(i).shape{phase}(j).weight = apertureInfoVect(updatedInfo.beam(i).shape{phase}(j).weightOffset)./updatedInfo.beam(i).shape{phase}(j).jacobiScale;
                
                if ~updatedInfo.propVMAT.continuousAperture
                    % extract left and right leaf positions from shape vector
                    vectorIx_L = updatedInfo.beam(i).shape{phase}(j).vectorOffset(1) + ((1:n)-1);
                    vectorIx_R = vectorIx_L+updatedInfo.totalNumOfLeafPairs*updatedInfo.numPhases;
                    leftLeafPos  = apertureInfoVect(vectorIx_L);
                    rightLeafPos = apertureInfoVect(vectorIx_R);
                    
                    % update information in shape structure
                    updatedInfo.beam(i).shape{phase}(j).leftLeafPos  = leftLeafPos;
                    updatedInfo.beam(i).shape{phase}(j).leftLeafPos_I = leftLeafPos;
                    updatedInfo.beam(i).shape{phase}(j).leftLeafPos_F = leftLeafPos;
                    updatedInfo.beam(i).shape{phase}(j).rightLeafPos = rightLeafPos;
                    updatedInfo.beam(i).shape{phase}(j).rightLeafPos_I = rightLeafPos;
                    updatedInfo.beam(i).shape{phase}(j).rightLeafPos_F = rightLeafPos;
                end
                
            else
                % this is an interpolated beam
                
                %MURate is interpolated between MURates of optimized apertures
                if phase == 1
                    updatedInfo.beam(i).gantryRot = 1./(updatedInfo.propVMAT.beam(i).fracFromLastDAO_gantryRot./updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).gantryRot+updatedInfo.propVMAT.beam(i).fracFromNextDAO_gantryRot./updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).gantryRot);
                    updatedInfo.beam(i).time = updatedInfo.propVMAT.beam(i).fluAngleBordersDiff./updatedInfo.beam(i).gantryRot;
                end
                updatedInfo.beam(i).shape{phase}(j).MURate = updatedInfo.propVMAT.beam(i).fracFromLastDAO_MU*updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(j).MURate+updatedInfo.propVMAT.beam(i).fracFromNextDAO_MU*updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(j).MURate;
                
                % calculate MU, weight
                updatedInfo.beam(i).shape{phase}(j).MU = updatedInfo.beam(i).shape{phase}(j).MURate.*updatedInfo.beam(i).time;
                updatedInfo.beam(i).shape{phase}(j).weight = updatedInfo.beam(i).shape{phase}(j).MU./updatedInfo.weightToMU;
                
                if ~updatedInfo.propVMAT.continuousAperture
                    
                    % obtain leaf positions at last DAO beam
                    vectorIx_lastDAOI = updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(j).vectorOffset + ((1:n)-1);
                    vectorIx_RF_last = vectorIx_lastDAOI+updatedInfo.totalNumOfLeafPairs*updatedInfo.numPhases;
                    leftLeafPos_last = apertureInfoVect(vectorIx_lastDAOI);
                    rightLeafPos_last = apertureInfoVect(vectorIx_RF_last);
                    
                    % obtain leaf positions at next DAO beam
                    vectorIx_LI_next = updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(j).vectorOffset + ((1:n)-1);
                    vectorIx_RI_next = vectorIx_LI_next+updatedInfo.totalNumOfLeafPairs*updatedInfo.numPhases;
                    leftLeafPos_next = apertureInfoVect(vectorIx_LI_next);
                    rightLeafPos_next = apertureInfoVect(vectorIx_RI_next);
                    
                    % interpolate leaf positions
                    leftLeafPos = updatedInfo.propVMAT.beam(i).fracFromLastDAO*leftLeafPos_last+(1-updatedInfo.propVMAT.beam(i).fracFromLastDAO)*leftLeafPos_next;
                    rightLeafPos = updatedInfo.propVMAT.beam(i).fracFromLastDAO*rightLeafPos_last+(1-updatedInfo.propVMAT.beam(i).fracFromLastDAO)*rightLeafPos_next;
                    
                    % update information in shape structure
                    updatedInfo.beam(i).shape{phase}(j).leftLeafPos  = leftLeafPos;
                    updatedInfo.beam(i).shape{phase}(j).leftLeafPos_I = leftLeafPos;
                    updatedInfo.beam(i).shape{phase}(j).leftLeafPos_F = leftLeafPos;
                    updatedInfo.beam(i).shape{phase}(j).rightLeafPos = rightLeafPos;
                    updatedInfo.beam(i).shape{phase}(j).rightLeafPos_I = rightLeafPos;
                    updatedInfo.beam(i).shape{phase}(j).rightLeafPos_F = rightLeafPos;
                end
            end
            
            if updatedInfo.runVMAT && updatedInfo.propVMAT.continuousAperture
                % leaf positions in the optimization vector are always
                % defined at the DAO angle boundaries
                % however, we must interpolate them at the fluence
                % calculation angle boundaries to calc fluence/dose
                
                % extract interpolation factors relating the optimized leaf
                % positions and the positions at the fluence boundaries
                fracFromLastDAOI_leafI = updatedInfo.propVMAT.beam(i).fracFromLastDAOI_leafI;
                fracFromLastDAOF_leafI = updatedInfo.propVMAT.beam(i).fracFromLastDAOF_leafI;
                fracFromNextDAOI_leafF = updatedInfo.propVMAT.beam(i).fracFromNextDAOI_leafF;
                fracFromNextDAOF_leafF = updatedInfo.propVMAT.beam(i).fracFromNextDAOF_leafF;
                
                % obtain initial and final leaf positions at last DAO beam
                vectorIx_lastDAOI       = updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(j).vectorOffset(1) + ((1:n)-1);
                vectorIx_lastDAOF       = updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(j).vectorOffset(2) + ((1:n)-1);
                leftLeafPos_lastDAOI    = apertureInfoVect(vectorIx_lastDAOI);
                rightLeafPos_lastDAOI   = apertureInfoVect(vectorIx_lastDAOI+updatedInfo.totalNumOfLeafPairs*updatedInfo.numPhases);
                leftLeafPos_lastDAOF    = apertureInfoVect(vectorIx_lastDAOF);
                rightLeafPos_lastDAOF   = apertureInfoVect(vectorIx_lastDAOF+updatedInfo.totalNumOfLeafPairs*updatedInfo.numPhases);
                
                % obtain initial and final leaf positions at next DAO beam
                vectorIx_nextDAOI       = updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(j).vectorOffset(1) + ((1:n)-1);
                vectorIx_nextDAOF       = updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(j).vectorOffset(2) + ((1:n)-1);
                leftLeafPos_nextDAOI    = apertureInfoVect(vectorIx_nextDAOI);
                rightLeafPos_nextDAOI   = apertureInfoVect(vectorIx_nextDAOI+updatedInfo.totalNumOfLeafPairs*updatedInfo.numPhases);
                leftLeafPos_nextDAOF    = apertureInfoVect(vectorIx_nextDAOF);
                rightLeafPos_nextDAOF   = apertureInfoVect(vectorIx_nextDAOF+updatedInfo.totalNumOfLeafPairs*updatedInfo.numPhases);
                
                % interpolate leaf positions
                updatedInfo.beam(i).shape{phase}(j).leftLeafPos_I   = fracFromLastDAOI_leafI.*leftLeafPos_lastDAOI+fracFromLastDAOF_leafI.*leftLeafPos_lastDAOF;
                updatedInfo.beam(i).shape{phase}(j).rightLeafPos_I  = fracFromLastDAOI_leafI.*rightLeafPos_lastDAOI+fracFromLastDAOF_leafI.*rightLeafPos_lastDAOF;
                
                updatedInfo.beam(i).shape{phase}(j).leftLeafPos_F   = fracFromNextDAOI_leafF.*leftLeafPos_nextDAOI+fracFromNextDAOF_leafF.*leftLeafPos_nextDAOF;
                updatedInfo.beam(i).shape{phase}(j).rightLeafPos_F  = fracFromNextDAOI_leafF.*rightLeafPos_nextDAOI+fracFromNextDAOF_leafF.*rightLeafPos_nextDAOF;
                
                % also update DAO leaf positions in struct
                if updatedInfo.propVMAT.beam(i).DAOBeam
                    if apertureInfo.propVMAT.beam(i).firstDAO
                        updatedInfo.beam(i).shape{phase}(j).leftLeafPos_I_DAO   = leftLeafPos_lastDAOI;
                        updatedInfo.beam(i).shape{phase}(j).rightLeafPos_I_DAO  = rightLeafPos_lastDAOI;
                    end
                    updatedInfo.beam(i).shape{phase}(j).leftLeafPos_F_DAO   = leftLeafPos_lastDAOF;
                    updatedInfo.beam(i).shape{phase}(j).rightLeafPos_F_DAO  = rightLeafPos_lastDAOF;
                end
            end
        end
    end
    
    % calculate all bixel weights and gradients for this gantry angle
    % ONLY PASS THE ARRAYS THAT WE ACTUALLY NEED FOR THIS ANGLE
    % REDUCE THE SIZE OF ARRAYS THAT ARE ACCESSED IN
    % matRad_bixWeightAndGrad
    
    % weights
    bixWeightAndGradBase_lastDose_angle.w(:) = {zeros(apertureInfo.beam(i).lastNumBixels,1)};
    bixWeightAndGradBase_nextDose_angle.w(:) = {zeros(apertureInfo.beam(i).nextNumBixels,1)};
    % jacobian
    bixWeightAndGradBase_lastDose_angle.bixelJApVec_vec(:) = {zeros(updatedInfo.beam(i).bixelJApVecLastDose_sz,1)};
    bixWeightAndGradBase_nextDose_angle.bixelJApVec_vec(:) = {zeros(updatedInfo.beam(i).bixelJApVecNextDose_sz,1)};
    % vector indices
    bixWeightAndGradBase_lastDose_angle.bixelJApVec_i(:) = {zeros(updatedInfo.beam(i).bixelJApVecLastDose_sz,1,iBitStr)};
    bixWeightAndGradBase_nextDose_angle.bixelJApVec_i(:) = {zeros(updatedInfo.beam(i).bixelJApVecNextDose_sz,1,iBitStr)};
    % bixel indices
    bixWeightAndGradBase_lastDose_angle.bixelJApVec_j(:) = {zeros(updatedInfo.beam(i).bixelJApVecLastDose_sz,1,jBitStr)};
    bixWeightAndGradBase_nextDose_angle.bixelJApVec_j(:) = {zeros(updatedInfo.beam(i).bixelJApVecNextDose_sz,1,jBitStr)};
    % offset
    bixWeightAndGradBase_lastDose_angle.bixelJApVec_offset(:) = {0};
    bixWeightAndGradBase_nextDose_angle.bixelJApVec_offset(:) = {0};
    
    % make arcI and arcF structs
    bixWeightAndGrad_angle.arcI.lastDose = bixWeightAndGradBase_lastDose_angle;
    bixWeightAndGrad_angle.arcI.nextDose = bixWeightAndGradBase_nextDose_angle;
    bixWeightAndGrad_angle.arcF.lastDose = bixWeightAndGradBase_lastDose_angle;
    bixWeightAndGrad_angle.arcF.nextDose = bixWeightAndGradBase_nextDose_angle;
    
    % initialize bixel gradients
    updatedInfo.bixelJApVec((i-1).*updatedInfo.numPhases+(1:updatedInfo.numPhases)) = {sparse(updatedInfo.beam(i).numUniqueVar,updatedInfo.beam(i).lastNumBixels+updatedInfo.beam(i).nextNumBixels)};
    
    % restrict the bixelJApVec jacobians to only the variables that can
    % possibly have an effect on the weights of that gantry angle
    [updatedInfo,bixWeightAndGrad_angle] = matRad_bixWeightAndGradWrapper(updatedInfo,i,bixWeightAndGrad_angle);
    
    % get bixel and fixel indices for last/next dose beams
    lastBixelIndMap = apertureInfo.beam(i).lastBixelIndMap(~isnan(apertureInfo.beam(i).lastBixelIndMap));
    nextBixelIndMap = apertureInfo.beam(i).nextBixelIndMap(~isnan(apertureInfo.beam(i).nextBixelIndMap));
    lastFixelIndMap = apertureInfo.beam(i).lastFixelIndMap(~isnan(apertureInfo.beam(i).lastFixelIndMap));
    nextFixelIndMap = apertureInfo.beam(i).nextFixelIndMap(~isnan(apertureInfo.beam(i).nextFixelIndMap));
    lastBixelIndMap_JApVec  = lastBixelIndMap-max(lastBixelIndMap(:))+apertureInfo.beam(i).lastNumBixels;
    nextBixelIndMap_JApVec  = nextBixelIndMap-max(nextBixelIndMap(:))+apertureInfo.beam(i).lastNumBixels+apertureInfo.beam(i).nextNumBixels;
    
    % loop over initial and final phases to put weights and gradients in
    % their proper place
    for phase_I = 1:updatedInfo.numPhases
        
        for phase_F = 1:updatedInfo.numPhases
            % DO BOTH LOOPS, PHASE_I AND PHASE_F, OUT HERE
            % have arcI.w be of length numPhases*numPhases
            % let it contain not the sum of all final phases but each final
            % phase individually
            % THEN we can use this for the variance!!
            
            cellInd = (phase_I-1).*updatedInfo.numPhases+phase_F;
            
            % split fixel weights into arcI/arcF and last/next dose beam
            updatedInfo.arcI.lastDose.fixelWeights{cellInd}(lastFixelIndMap) = bixWeightAndGrad_angle.arcI.lastDose.w{cellInd};
            updatedInfo.arcI.nextDose.fixelWeights{cellInd}(nextFixelIndMap) = bixWeightAndGrad_angle.arcI.nextDose.w{cellInd};
            updatedInfo.arcF.lastDose.fixelWeights{cellInd}(lastFixelIndMap) = bixWeightAndGrad_angle.arcF.lastDose.w{cellInd};
            updatedInfo.arcF.nextDose.fixelWeights{cellInd}(nextFixelIndMap) = bixWeightAndGrad_angle.arcF.nextDose.w{cellInd};
            
            % sum real bixel weights (both initial and final phase)
            updatedInfo.bixelWeights{phase_I}(lastBixelIndMap) = updatedInfo.bixelWeights{phase_I}(lastBixelIndMap)+bixWeightAndGrad_angle.arcI.lastDose.w{cellInd};
            updatedInfo.bixelWeights{phase_I}(nextBixelIndMap) = updatedInfo.bixelWeights{phase_I}(nextBixelIndMap)+bixWeightAndGrad_angle.arcI.nextDose.w{cellInd};
            updatedInfo.bixelWeights{phase_F}(lastBixelIndMap) = updatedInfo.bixelWeights{phase_F}(lastBixelIndMap)+bixWeightAndGrad_angle.arcF.lastDose.w{cellInd};
            updatedInfo.bixelWeights{phase_F}(nextBixelIndMap) = updatedInfo.bixelWeights{phase_F}(nextBixelIndMap)+bixWeightAndGrad_angle.arcF.nextDose.w{cellInd};
            
            % split fixel gradients into arcI/arcF and last/next dose beam
            updatedInfo.arcI.lastDose.bixelJApVec{(i-1).*updatedInfo.numPhases.^2+cellInd} = sparse(double(bixWeightAndGrad_angle.arcI.lastDose.bixelJApVec_i{cellInd}),double(bixWeightAndGrad_angle.arcI.lastDose.bixelJApVec_j{cellInd}),bixWeightAndGrad_angle.arcI.lastDose.bixelJApVec_vec{cellInd},updatedInfo.beam(i).numUniqueVar,updatedInfo.beam(i).lastNumBixels);
            updatedInfo.arcI.nextDose.bixelJApVec{(i-1).*updatedInfo.numPhases.^2+cellInd} = sparse(double(bixWeightAndGrad_angle.arcI.nextDose.bixelJApVec_i{cellInd}),double(bixWeightAndGrad_angle.arcI.nextDose.bixelJApVec_j{cellInd}),bixWeightAndGrad_angle.arcI.nextDose.bixelJApVec_vec{cellInd},updatedInfo.beam(i).numUniqueVar,updatedInfo.beam(i).nextNumBixels);
            updatedInfo.arcF.lastDose.bixelJApVec{(i-1).*updatedInfo.numPhases.^2+cellInd} = sparse(double(bixWeightAndGrad_angle.arcF.lastDose.bixelJApVec_i{cellInd}),double(bixWeightAndGrad_angle.arcF.lastDose.bixelJApVec_j{cellInd}),bixWeightAndGrad_angle.arcF.lastDose.bixelJApVec_vec{cellInd},updatedInfo.beam(i).numUniqueVar,updatedInfo.beam(i).lastNumBixels);
            updatedInfo.arcF.nextDose.bixelJApVec{(i-1).*updatedInfo.numPhases.^2+cellInd} = sparse(double(bixWeightAndGrad_angle.arcF.nextDose.bixelJApVec_i{cellInd}),double(bixWeightAndGrad_angle.arcF.nextDose.bixelJApVec_j{cellInd}),bixWeightAndGrad_angle.arcF.nextDose.bixelJApVec_vec{cellInd},updatedInfo.beam(i).numUniqueVar,updatedInfo.beam(i).nextNumBixels);
            
            % sum real bixel gradients (both initial and final phase)
            updatedInfo.bixelJApVec{(i-1).*updatedInfo.numPhases+phase_I}(:,lastBixelIndMap_JApVec) = updatedInfo.bixelJApVec{(i-1).*updatedInfo.numPhases+phase_I}(:,lastBixelIndMap_JApVec)+updatedInfo.arcI.lastDose.bixelJApVec{(i-1).*updatedInfo.numPhases.^2+cellInd};
            updatedInfo.bixelJApVec{(i-1).*updatedInfo.numPhases+phase_I}(:,nextBixelIndMap_JApVec) = updatedInfo.bixelJApVec{(i-1).*updatedInfo.numPhases+phase_I}(:,nextBixelIndMap_JApVec)+updatedInfo.arcI.nextDose.bixelJApVec{(i-1).*updatedInfo.numPhases.^2+cellInd};
            updatedInfo.bixelJApVec{(i-1).*updatedInfo.numPhases+phase_F}(:,lastBixelIndMap_JApVec) = updatedInfo.bixelJApVec{(i-1).*updatedInfo.numPhases+phase_F}(:,lastBixelIndMap_JApVec)+updatedInfo.arcF.lastDose.bixelJApVec{(i-1).*updatedInfo.numPhases.^2+cellInd};
            updatedInfo.bixelJApVec{(i-1).*updatedInfo.numPhases+phase_F}(:,nextBixelIndMap_JApVec) = updatedInfo.bixelJApVec{(i-1).*updatedInfo.numPhases+phase_F}(:,nextBixelIndMap_JApVec)+updatedInfo.arcF.nextDose.bixelJApVec{(i-1).*updatedInfo.numPhases.^2+cellInd};
            
        end
    end
end

end