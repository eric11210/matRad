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

bixWeightAndGradBase.w = cell(updatedInfo.numPhases,1);
bixWeightAndGradBase.w(:) = {zeros(updatedInfo.totalNumOfBixels,1)};

if updatedInfo.runVMAT && ~all([updatedInfo.propVMAT.beam.DAOBeam])
    j = 1;
    for phase = 1:updatedInfo.numPhases
        for i = 1:numel(updatedInfo.beam)
            if updatedInfo.propVMAT.beam(i).DAOBeam
                % update the shape weight
                % rescale the weight from the vector using the previous
                % iteration scaling factor
                updatedInfo.beam(i).shape{phase}(j).weight = apertureInfoVect(updatedInfo.beam(i).shape{phase}(j).weightOffset)./updatedInfo.beam(i).shape{phase}(j).jacobiScale;
                
                updatedInfo.beam(i).shape{phase}(j).MU = updatedInfo.beam(i).shape{phase}(j).weight*updatedInfo.weightToMU;
                if phase == 1
                    updatedInfo.beam(i).time = apertureInfoVect((updatedInfo.totalNumOfShapes+updatedInfo.totalNumOfLeafPairs*2)*updatedInfo.numPhases+updatedInfo.propVMAT.beam(i).DAOIndex)*updatedInfo.propVMAT.beam(i).timeFacCurr;
                    updatedInfo.beam(i).gantryRot = updatedInfo.propVMAT.beam(i).doseAngleBordersDiff/updatedInfo.beam(i).time;
                end
                updatedInfo.beam(i).shape{phase}(j).MURate = updatedInfo.beam(i).shape{phase}(j).MU./updatedInfo.beam(i).time;
            end
        end
    end
end


% Jacobian matrix to be used in the DAO gradient function
% this tells us the gradient of a particular bixel with respect to an
% element in the apertureVector (aperture weight or leaf position)
% store as a vector for now, convert to sparse matrix later

optBixelFactor = 7;
% For optimized beams: 7 = (1 from weights) + (3 from left leaf positions (I, M, and F)) + (3 from
% right leaf positions (I, M, and F))

if updatedInfo.runVMAT
    intBixelFactor = 2*optBixelFactor+2;
    % For interpolated beams: multiply this number times 2 (influenced by the
    % one before and the one after), then add 2 (influenced by the time of the
    % times before and after)
else
    intBixelFactor = 2*optBixelFactor;
    % For interpolated beams: multiply this number times 2 (influenced by the
    % one before and the one after)
end

% for the time (probability) gradients
optBixelFactor = optBixelFactor+updatedInfo.totalNumOfShapes;
intBixelFactor = intBixelFactor+updatedInfo.totalNumOfShapes;
%bixelJApVec_sz = (updatedInfo.totalNumOfOptBixels*optBixelFactor+(updatedInfo.totalNumOfBixels-updatedInfo.totalNumOfOptBixels)*intBixelFactor)*updatedInfo.numPhases*2;
bixelJApVec_sz = (updatedInfo.totalNumOfOptBixels*optBixelFactor+(updatedInfo.totalNumOfBixels-updatedInfo.totalNumOfOptBixels)*intBixelFactor)*updatedInfo.numPhases;

% jacobian
bixWeightAndGradBase.bixelJApVec_vec = cell(updatedInfo.numPhases,1);
bixWeightAndGradBase.bixelJApVec_vec(:) = {zeros(1,bixelJApVec_sz)};
% vector indices
bixWeightAndGradBase.bixelJApVec_i = cell(updatedInfo.numPhases,1);
bixWeightAndGradBase.bixelJApVec_i(:) = {nan(1,bixelJApVec_sz)};
% bixel indices
bixWeightAndGradBase.bixelJApVec_j = cell(updatedInfo.numPhases,1);
bixWeightAndGradBase.bixelJApVec_j(:) = {zeros(1,bixelJApVec_sz)};
% offset
bixWeightAndGradBase.bixelJApVec_offset = cell(updatedInfo.numPhases,1);
bixWeightAndGradBase.bixelJApVec_offset(:) = {0};

% make arcI and arcF structs
bixWeightAndGrad.arcI = bixWeightAndGradBase;
bixWeightAndGrad.arcF = bixWeightAndGradBase;

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
                
                if updatedInfo.runVMAT
                    updatedInfo.beam(i).shape{phase}(j).MU = updatedInfo.beam(i).shape{phase}(j).weight*updatedInfo.weightToMU;
                    if phase == 1
                        updatedInfo.beam(i).time = apertureInfoVect((updatedInfo.totalNumOfShapes+updatedInfo.totalNumOfLeafPairs*2)*updatedInfo.numPhases+updatedInfo.propVMAT.beam(i).DAOIndex)*updatedInfo.propVMAT.beam(i).timeFacCurr;
                        updatedInfo.beam(i).gantryRot = updatedInfo.propVMAT.beam(i).doseAngleBordersDiff/updatedInfo.beam(i).time;
                    end
                    updatedInfo.beam(i).shape{phase}(j).MURate = updatedInfo.beam(i).shape{phase}(j).MU./updatedInfo.beam(i).time;
                end
                
                if ~updatedInfo.propVMAT.continuousAperture
                    % extract left and right leaf positions from shape vector
                    vectorIx_L = updatedInfo.beam(i).shape{phase}(j).vectorOffset + ((1:n)-1);
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
                else
                    % extract left and right leaf positions from shape vector
                    vectorIx_LI = updatedInfo.beam(i).shape{phase}(j).vectorOffset(1) + ((1:n)-1);
                    vectorIx_RI = vectorIx_LI+updatedInfo.totalNumOfLeafPairs*updatedInfo.numPhases;
                    leftLeafPos_I = apertureInfoVect(vectorIx_LI);
                    rightLeafPos_I = apertureInfoVect(vectorIx_RI);
                    
                    vectorIx_LF = updatedInfo.beam(i).shape{phase}(j).vectorOffset(2) + ((1:n)-1);
                    vectorIx_RF = vectorIx_LF+updatedInfo.totalNumOfLeafPairs*updatedInfo.numPhases;
                    leftLeafPos_F = apertureInfoVect(vectorIx_LF);
                    rightLeafPos_F = apertureInfoVect(vectorIx_RF);
                    
                    % update information in shape structure
                    updatedInfo.beam(i).shape{phase}(j).leftLeafPos_I  = leftLeafPos_I;
                    updatedInfo.beam(i).shape{phase}(j).rightLeafPos_I = rightLeafPos_I;
                    
                    updatedInfo.beam(i).shape{phase}(j).leftLeafPos_F  = leftLeafPos_F;
                    updatedInfo.beam(i).shape{phase}(j).rightLeafPos_F = rightLeafPos_F;
                end
                
            else
                % this is an interpolated beam
                
                %MURate is interpolated between MURates of optimized apertures
                if phase == 1
                    updatedInfo.beam(i).gantryRot = 1./(updatedInfo.propVMAT.beam(i).timeFracFromLastDAO./updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).gantryRot+updatedInfo.propVMAT.beam(i).timeFracFromNextDAO./updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).gantryRot);
                    updatedInfo.beam(i).time = updatedInfo.propVMAT.beam(i).doseAngleBordersDiff./updatedInfo.beam(i).gantryRot;
                end
                updatedInfo.beam(i).shape{phase}(j).MURate = updatedInfo.propVMAT.beam(i).fracFromLastDAO*updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(j).MURate+(1-updatedInfo.propVMAT.beam(i).fracFromLastDAO)*updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(j).MURate;
                
                % calculate MU, weight
                updatedInfo.beam(i).shape{phase}(j).MU = updatedInfo.beam(i).shape{phase}(j).MURate.*updatedInfo.beam(i).time;
                updatedInfo.beam(i).shape{phase}(j).weight = updatedInfo.beam(i).shape{phase}(j).MU./updatedInfo.weightToMU;
                
                if ~updatedInfo.propVMAT.continuousAperture
                    
                    % obtain leaf positions at last DAO beam
                    vectorIx_LF_last = updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(j).vectorOffset + ((1:n)-1);
                    vectorIx_RF_last = vectorIx_LF_last+updatedInfo.totalNumOfLeafPairs*updatedInfo.numPhases;
                    leftLeafPos_last = apertureInfoVect(vectorIx_LF_last);
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
                else
                    
                    fracFromLastOptI = updatedInfo.propVMAT.beam(i).fracFromLastDAO_I*ones(n,1);
                    fracFromLastOptF = updatedInfo.propVMAT.beam(i).fracFromLastDAO_F*ones(n,1);
                    fracFromNextOptI = updatedInfo.propVMAT.beam(i).fracFromNextDAO_I*ones(n,1);
                    fracFromNextOptF = updatedInfo.propVMAT.beam(i).fracFromNextDAO_F*ones(n,1);
                    
                    % obtain leaf positions at last DAO beam
                    vectorIx_LF_last = updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(j).vectorOffset(2) + ((1:n)-1);
                    vectorIx_RF_last = vectorIx_LF_last+updatedInfo.totalNumOfLeafPairs*updatedInfo.numPhases;
                    leftLeafPos_F_last = apertureInfoVect(vectorIx_LF_last);
                    rightLeafPos_F_last = apertureInfoVect(vectorIx_RF_last);
                    
                    % obtain leaf positions at next DAO beam
                    vectorIx_LI_next = updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(j).vectorOffset(1) + ((1:n)-1);
                    vectorIx_RI_next = vectorIx_LI_next+updatedInfo.totalNumOfLeafPairs*updatedInfo.numPhases;
                    leftLeafPos_I_next = apertureInfoVect(vectorIx_LI_next);
                    rightLeafPos_I_next = apertureInfoVect(vectorIx_RI_next);
                    
                    % interpolate leaf positions
                    updatedInfo.beam(i).shape{phase}(j).leftLeafPos_I = fracFromLastOptI.*leftLeafPos_F_last+fracFromNextOptI.*leftLeafPos_I_next;
                    updatedInfo.beam(i).shape{phase}(j).rightLeafPos_I = fracFromLastOptI.*rightLeafPos_F_last+fracFromNextOptI.*rightLeafPos_I_next;
                    
                    updatedInfo.beam(i).shape{phase}(j).leftLeafPos_F = fracFromLastOptF.*leftLeafPos_F_last+fracFromNextOptF.*leftLeafPos_I_next;
                    updatedInfo.beam(i).shape{phase}(j).rightLeafPos_F = fracFromLastOptF.*rightLeafPos_F_last+fracFromNextOptF.*rightLeafPos_I_next;
                end
            end
        end
    end
    
    % calculate all bixel weights and gradients for this gantry angle
    [updatedInfo,bixWeightAndGrad] = matRad_bixWeightAndGradWrapper(updatedInfo,i,bixWeightAndGrad);
end

%% save bixelWeight, apertureVector, and Jacobian between the two
updatedInfo.arcI.bixelWeights = bixWeightAndGrad.arcI.w;
updatedInfo.arcF.bixelWeights = bixWeightAndGrad.arcF.w;
updatedInfo.bixelWeights = cell(updatedInfo.numPhases,1);

updatedInfo.apertureVector = apertureInfoVect;

updatedInfo.arcI.bixelJApVec = cell(updatedInfo.numPhases,1);
updatedInfo.arcF.bixelJApVec = cell(updatedInfo.numPhases,1);
updatedInfo.bixelJApVec = cell(updatedInfo.numPhases,1);
for phase = 1:updatedInfo.numPhases
    
    deleteInd = isnan(bixWeightAndGrad.arcI.bixelJApVec_i{phase});
    bixWeightAndGrad.arcI.bixelJApVec_i{phase}(deleteInd) = [];
    bixWeightAndGrad.arcI.bixelJApVec_j{phase}(deleteInd) = [];
    bixWeightAndGrad.arcI.bixelJApVec_vec{phase}(deleteInd) = [];
    
    deleteInd = isnan(bixWeightAndGrad.arcF.bixelJApVec_i{phase});
    bixWeightAndGrad.arcF.bixelJApVec_i{phase}(deleteInd) = [];
    bixWeightAndGrad.arcF.bixelJApVec_j{phase}(deleteInd) = [];
    bixWeightAndGrad.arcF.bixelJApVec_vec{phase}(deleteInd) = [];
    
    updatedInfo.arcI.bixelJApVec{phase} = sparse(bixWeightAndGrad.arcI.bixelJApVec_i{phase},bixWeightAndGrad.arcI.bixelJApVec_j{phase},bixWeightAndGrad.arcI.bixelJApVec_vec{phase},numel(apertureInfoVect),updatedInfo.totalNumOfBixels);
    updatedInfo.arcF.bixelJApVec{phase} = sparse(bixWeightAndGrad.arcF.bixelJApVec_i{phase},bixWeightAndGrad.arcF.bixelJApVec_j{phase},bixWeightAndGrad.arcF.bixelJApVec_vec{phase},numel(apertureInfoVect),updatedInfo.totalNumOfBixels);
    
    % sum both arcs
    updatedInfo.bixelWeights{phase}    = updatedInfo.arcI.bixelWeights{phase}+updatedInfo.arcF.bixelWeights{phase};
    updatedInfo.bixelJApVec{phase}     = updatedInfo.arcI.bixelJApVec{phase}+updatedInfo.arcF.bixelJApVec{phase};
end

end