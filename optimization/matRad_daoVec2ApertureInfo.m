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

shapeInd = 1;

indVect = NaN*ones(2*apertureInfo.doseTotalNumOfLeafPairs*apertureInfo.numPhases,1);
offset = 0;

% options for bixel and Jacobian calculation
mlcOptions.bixelWidth = apertureInfo.bixelWidth;
calcOptions.continuousAperture = updatedInfo.propVMAT.continuousAperture;
vectorIndices.totalNumOfShapes = apertureInfo.totalNumOfShapes;
vectorIndices.timeOffset = (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases;

w = cell(apertureInfo.numPhases,1);
w(:) = {zeros(apertureInfo.totalNumOfBixels,1)};

% helper function to cope with numerical instabilities through rounding
round2 = @(a,b) round(a*10^b)/10^b;

if updatedInfo.runVMAT && ~all([updatedInfo.propVMAT.beam.DAOBeam])
    j = 1;
    for phase = 1:apertureInfo.numPhases
        for i = 1:numel(updatedInfo.beam)
            if updatedInfo.propVMAT.beam(i).DAOBeam
                % update the shape weight
                % rescale the weight from the vector using the previous
                % iteration scaling factor
                updatedInfo.beam(i).shape{phase}(j).weight = apertureInfoVect(shapeInd)./updatedInfo.beam(i).shape{phase}(j).jacobiScale;
                
                updatedInfo.beam(i).shape{phase}(j).MU = updatedInfo.beam(i).shape{phase}(j).weight*updatedInfo.weightToMU;
                if phase == 1
                    updatedInfo.beam(i).time = apertureInfoVect((updatedInfo.totalNumOfShapes+updatedInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases+shapeInd)*updatedInfo.propVMAT.beam(i).timeFacCurr;
                    updatedInfo.beam(i).gantryRot = updatedInfo.propVMAT.beam(i).doseAngleBordersDiff/updatedInfo.beam(i).time;
                end
                updatedInfo.beam(i).shape{phase}(j).MURate = updatedInfo.beam(i).shape{phase}(j).MU./updatedInfo.beam(i).time;
                
                shapeInd = shapeInd+1;
            end
        end
    end
    shapeInd = 1;
end


% Jacobian matrix to be used in the DAO gradient function
% this tells us the gradient of a particular bixel with respect to an
% element in the apertureVector (aperture weight or leaf position)
% store as a vector for now, convert to sparse matrix later

if updatedInfo.propVMAT.continuousAperture
    optBixelFactor = 5;
    % For optimized beams: 5 = (1 from weights) + (2 from left leaf positions (I and F)) + (2 from
    % right leaf positions (I and F))
else
    optBixelFactor = 3;
    % For optimized beams: 3 = (1 from weights) + (1 from left leaf positions) + (1 from
    % right leaf positions )
end

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

bixelJApVec_sz = updatedInfo.totalNumOfOptBixels*optBixelFactor+(updatedInfo.totalNumOfBixels-updatedInfo.totalNumOfOptBixels)*intBixelFactor;

bixelJApVec_vec = cell(apertureInfo.numPhases,1);
bixelJApVec_vec(:) = {zeros(1,bixelJApVec_sz)};

% vector indices
bixelJApVec_i = cell(apertureInfo.numPhases,1);
bixelJApVec_i(:) = {zeros(1,bixelJApVec_sz)};
% bixel indices
bixelJApVec_j = cell(apertureInfo.numPhases,1);
bixelJApVec_j(:) = {zeros(1,bixelJApVec_sz)};
% offset
bixelJApVec_offset = 0;



%% update the shapeMaps
% here the new colimator positions are used to create new shapeMaps that
% now include decimal values instead of binary

% loop over all beams
for phase = 1:apertureInfo.numPhases
    
    for i = 1:numel(updatedInfo.beam)
        
        %posOfRightCornerPixel = apertureInfo.beam(i).posOfCornerBixel(1) + (size(apertureInfo.beam(i).bixelIndMap,2)-1)*apertureInfo.bixelWidth;
        
        % pre compute left and right bixel edges
        edges_l = updatedInfo.beam(i).posOfCornerBixel(1)...
            + ([1:size(apertureInfo.beam(i).bixelIndMap,2)]-1-1/2)*updatedInfo.bixelWidth;
        edges_r = updatedInfo.beam(i).posOfCornerBixel(1)...
            + ([1:size(apertureInfo.beam(i).bixelIndMap,2)]-1+1/2)*updatedInfo.bixelWidth;
        
        % get dimensions of 2d matrices that store shape/bixel information
        n = apertureInfo.beam(i).numOfActiveLeafPairs;
        
        % loop over all shapes
        if updatedInfo.runVMAT
            numOfShapes = 1;
        else
            numOfShapes = updatedInfo.beam(i).numOfShapes;
        end
        
        mlcOptions.lim_l = apertureInfo.beam(i).lim_l;
        mlcOptions.lim_r = apertureInfo.beam(i).lim_r;
        mlcOptions.edges_l = edges_l;
        mlcOptions.edges_r = edges_r;
        mlcOptions.centres = (edges_l+edges_r)/2;
        mlcOptions.widths = edges_r-edges_l;
        mlcOptions.n = n;
        mlcOptions.numBix = size(apertureInfo.beam(i).bixelIndMap,2);
        mlcOptions.bixelIndMap = apertureInfo.beam(i).bixelIndMap;
        calcOptions.DAOBeam = updatedInfo.propVMAT.beam(i).DAOBeam;
        
        for j = 1:numOfShapes
            
            if ~updatedInfo.runVMAT || updatedInfo.propVMAT.beam(i).DAOBeam
                % either this is not VMAT, or if it is VMAT, this is a DAO beam
                
                % update the shape weight
                updatedInfo.beam(i).shape{phase}(j).weight = apertureInfoVect(shapeInd)./updatedInfo.beam(i).shape{phase}(j).jacobiScale;
                
                if updatedInfo.runVMAT
                    updatedInfo.beam(i).shape{phase}(j).MU = updatedInfo.beam(i).shape{phase}(j).weight*updatedInfo.weightToMU;
                    if phase == 1
                        updatedInfo.beam(i).time = apertureInfoVect((updatedInfo.totalNumOfShapes+updatedInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases+shapeInd)*updatedInfo.propVMAT.beam(i).timeFacCurr;
                        updatedInfo.beam(i).gantryRot = updatedInfo.propVMAT.beam(i).doseAngleBordersDiff/updatedInfo.beam(i).time;
                    end
                    updatedInfo.beam(i).shape{phase}(j).MURate = updatedInfo.beam(i).shape{phase}(j).MU./updatedInfo.beam(i).time;
                end
                
                if ~updatedInfo.propVMAT.continuousAperture
                    % extract left and right leaf positions from shape vector
                    vectorIx_L = updatedInfo.beam(i).shape{phase}(j).vectorOffset + ([1:n]-1);
                    vectorIx_LI = vectorIx_L;
                    vectorIx_LF = vectorIx_L;
                    vectorIx_R = vectorIx_L+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
                    vectorIx_RI = vectorIx_R;
                    vectorIx_RF = vectorIx_R;
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
                    vectorIx_LI = updatedInfo.beam(i).shape{phase}(j).vectorOffset(1) + ([1:n]-1);
                    vectorIx_RI = vectorIx_LI+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
                    leftLeafPos_I = apertureInfoVect(vectorIx_LI);
                    rightLeafPos_I = apertureInfoVect(vectorIx_RI);
                    
                    vectorIx_LF = updatedInfo.beam(i).shape{phase}(j).vectorOffset(2) + ([1:n]-1);
                    vectorIx_RF = vectorIx_LF+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
                    leftLeafPos_F = apertureInfoVect(vectorIx_LF);
                    rightLeafPos_F = apertureInfoVect(vectorIx_RF);
                    
                    % update information in shape structure
                    updatedInfo.beam(i).shape{phase}(j).leftLeafPos_I  = leftLeafPos_I;
                    updatedInfo.beam(i).shape{phase}(j).rightLeafPos_I = rightLeafPos_I;
                    
                    updatedInfo.beam(i).shape{phase}(j).leftLeafPos_F  = leftLeafPos_F;
                    updatedInfo.beam(i).shape{phase}(j).rightLeafPos_F = rightLeafPos_F;
                end
                
                % increment shape index
                shapeInd = shapeInd +1;
            else
                % this is an interpolated beam
                
                %MURate is interpolated between MURates of optimized apertures
                updatedInfo.beam(i).shape{phase}(j).MURate = updatedInfo.propVMAT.beam(i).fracFromLastDAO*updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(j).MURate+(1-updatedInfo.propVMAT.beam(i).fracFromLastDAO)*updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(j).MURate;
                updatedInfo.beam(i).gantryRot = 1./(updatedInfo.propVMAT.beam(i).timeFracFromLastDAO./updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).gantryRot+updatedInfo.propVMAT.beam(i).timeFracFromNextDAO./updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).gantryRot);
                updatedInfo.beam(i).time = updatedInfo.propVMAT.beam(i).doseAngleBordersDiff./updatedInfo.beam(i).gantryRot;
                
                % calculate MU, weight
                updatedInfo.beam(i).shape{phase}(j).MU = updatedInfo.beam(i).shape{phase}(j).MURate.*updatedInfo.beam(i).time;
                updatedInfo.beam(i).shape{phase}(j).weight = updatedInfo.beam(i).shape{phase}(j).MU./updatedInfo.weightToMU;
                
                if ~updatedInfo.propVMAT.continuousAperture
                    
                    fracFromLastOpt = updatedInfo.propVMAT.beam(i).fracFromLastDAO;
                    fracFromLastOptI = updatedInfo.propVMAT.beam(i).fracFromLastDAO*ones(n,1);
                    fracFromLastOptF = updatedInfo.propVMAT.beam(i).fracFromLastDAO*ones(n,1);
                    fracFromNextOptI = (1-updatedInfo.propVMAT.beam(i).fracFromLastDAO)*ones(n,1);
                    fracFromNextOptF = (1-updatedInfo.propVMAT.beam(i).fracFromLastDAO)*ones(n,1);
                    
                    % obtain leaf positions at last DAO beam
                    vectorIx_LF_last = updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(j).vectorOffset + ([1:n]-1);
                    vectorIx_RF_last = vectorIx_LF_last+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
                    leftLeafPos_last = apertureInfoVect(vectorIx_LF_last);
                    rightLeafPos_last = apertureInfoVect(vectorIx_RF_last);
                    
                    % obtain leaf positions at next DAO beam
                    vectorIx_LI_next = updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(j).vectorOffset + ([1:n]-1);
                    vectorIx_RI_next = vectorIx_LI_next+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
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
                    
                    fracFromLastOpt = updatedInfo.propVMAT.beam(i).fracFromLastDAO;
                    fracFromLastOptI = updatedInfo.propVMAT.beam(i).fracFromLastDAO_I*ones(n,1);
                    fracFromLastOptF = updatedInfo.propVMAT.beam(i).fracFromLastDAO_F*ones(n,1);
                    fracFromNextOptI = updatedInfo.propVMAT.beam(i).fracFromNextDAO_I*ones(n,1);
                    fracFromNextOptF = updatedInfo.propVMAT.beam(i).fracFromNextDAO_F*ones(n,1);
                    
                    % obtain leaf positions at last DAO beam
                    vectorIx_LF_last = updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(j).vectorOffset(2) + ([1:n]-1);
                    vectorIx_RF_last = vectorIx_LF_last+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
                    leftLeafPos_F_last = apertureInfoVect(vectorIx_LF_last);
                    rightLeafPos_F_last = apertureInfoVect(vectorIx_RF_last);
                    
                    % obtain leaf positions at next DAO beam
                    vectorIx_LI_next = updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(j).vectorOffset(1) + ([1:n]-1);
                    vectorIx_RI_next = vectorIx_LI_next+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
                    leftLeafPos_I_next = apertureInfoVect(vectorIx_LI_next);
                    rightLeafPos_I_next = apertureInfoVect(vectorIx_RI_next);
                    
                    % interpolate leaf positions
                    updatedInfo.beam(i).shape{phase}(j).leftLeafPos_I = fracFromLastOptI.*leftLeafPos_F_last+fracFromNextOptI.*leftLeafPos_I_next;
                    updatedInfo.beam(i).shape{phase}(j).rightLeafPos_I = fracFromLastOptI.*rightLeafPos_F_last+fracFromNextOptI.*rightLeafPos_I_next;
                    
                    updatedInfo.beam(i).shape{phase}(j).leftLeafPos_F = fracFromLastOptF.*leftLeafPos_F_last+fracFromNextOptF.*leftLeafPos_I_next;
                    updatedInfo.beam(i).shape{phase}(j).rightLeafPos_F = fracFromLastOptF.*rightLeafPos_F_last+fracFromNextOptF.*rightLeafPos_I_next;
                end
            end
            
            if ~updatedInfo.propVMAT.continuousAperture
                %% discrete aperture approximation
                
                % rounding for numerical stability
                leftLeafPos  = round2(leftLeafPos,10);
                rightLeafPos = round2(rightLeafPos,10);
                
                % check overshoot of leaf positions
                leftLeafPos(leftLeafPos <= apertureInfo.beam(i).lim_l) = apertureInfo.beam(i).lim_l(leftLeafPos <= apertureInfo.beam(i).lim_l);
                rightLeafPos(rightLeafPos <= apertureInfo.beam(i).lim_l) = apertureInfo.beam(i).lim_l(rightLeafPos <= apertureInfo.beam(i).lim_l);
                leftLeafPos(leftLeafPos >= apertureInfo.beam(i).lim_r) = apertureInfo.beam(i).lim_r(leftLeafPos >= apertureInfo.beam(i).lim_r);
                rightLeafPos(rightLeafPos >= apertureInfo.beam(i).lim_r) = apertureInfo.beam(i).lim_r(rightLeafPos >= apertureInfo.beam(i).lim_r);
                
                xPosIndLeftLeaf  = round((leftLeafPos - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth + 1);
                xPosIndRightLeaf = round((rightLeafPos - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth + 1);
                
                %
                xPosIndLeftLeaf_lim  = floor((apertureInfo.beam(i).lim_l - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth+1);
                xPosIndRightLeaf_lim = ceil((apertureInfo.beam(i).lim_r - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth + 1);
                
                xPosIndLeftLeaf(xPosIndLeftLeaf <= xPosIndLeftLeaf_lim) = xPosIndLeftLeaf_lim(xPosIndLeftLeaf <= xPosIndLeftLeaf_lim)+1;
                xPosIndRightLeaf(xPosIndRightLeaf >= xPosIndRightLeaf_lim) = xPosIndRightLeaf_lim(xPosIndRightLeaf >= xPosIndRightLeaf_lim)-1;
                
                % check limits because of rounding off issues at maximum, i.e.,
                % enforce round(X.5) -> X
                % LeafPos can occasionally go slightly beyond lim_r, so changed
                % == check to >=
                xPosIndLeftLeaf(leftLeafPos >= apertureInfo.beam(i).lim_r) = round(...
                    .5 + (leftLeafPos(leftLeafPos >= apertureInfo.beam(i).lim_r) ...
                    - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth);
                
                xPosIndRightLeaf(rightLeafPos >= apertureInfo.beam(i).lim_r) = round(...
                    .5 + (rightLeafPos(rightLeafPos >= apertureInfo.beam(i).lim_r) ...
                    - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth);
                
                % find the bixel index that the leaves currently touch
                bixelIndLeftLeaf  = apertureInfo.beam(i).bixelIndMap((xPosIndLeftLeaf-1)*n+[1:n]');
                bixelIndRightLeaf = apertureInfo.beam(i).bixelIndMap((xPosIndRightLeaf-1)*n+[1:n]');
                
                if any(isnan(bixelIndLeftLeaf)) || any(isnan(bixelIndRightLeaf))
                    error('cannot map leaf position to bixel index');
                end
                
                % store information in index vector for gradient calculation
                indVect(offset+(1:n)) = bixelIndLeftLeaf;
                indVect(offset+apertureInfo.doseTotalNumOfLeafPairs+(1:n)) = bixelIndRightLeaf;
                offset = offset+n;
                
                % calculate opening fraction for every bixel in shape to construct
                % bixel weight vector
                
                coveredByLeftLeaf  = bsxfun(@minus,leftLeafPos,edges_l)  / updatedInfo.bixelWidth;
                coveredByRightLeaf = bsxfun(@minus,edges_r,rightLeafPos) / updatedInfo.bixelWidth;
                
                tempMap = 1 - (coveredByLeftLeaf  + abs(coveredByLeftLeaf))  / 2 ...
                    - (coveredByRightLeaf + abs(coveredByRightLeaf)) / 2;
                
                % find open bixels
                tempMapIx = tempMap > 0;
                
                currBixelIx = apertureInfo.beam(i).bixelIndMap(tempMapIx);
                %w{phase}(currBixelIx) = w{phase}(currBixelIx) + tempMap(tempMapIx)*updatedInfo.beam(i).shape{phase}(j).weight;
                
                % save the tempMap (we need to apply a positivity operator !)
                updatedInfo.beam(i).shape{phase}(j).shapeMap = (tempMap  + abs(tempMap))  / 2;
                
            end
            
            %% enter in variables and options
            
            variables.phase = phase;
            variables.weight = updatedInfo.beam(i).shape{phase}(j).weight;
            
            % is it necessary to distinguish the apertureMotion for the
            % calculation part? (it definitely is for the bixelJApVec)
            
            if updatedInfo.propVMAT.beam(i).DAOBeam
                % the weight here is going to depend on the beginning and
                % ending phases
                % change phase to phase_I and phase_F?
                
                variables.jacobiScale = apertureInfo.beam(i).shape{phase}(1).jacobiScale;
                variables.leftLeafPos_I = updatedInfo.beam(i).shape{phase}(j).leftLeafPos_I;
                variables.leftLeafPos_F = updatedInfo.beam(i).shape{phase}(j).leftLeafPos_F;
                variables.rightLeafPos_I = updatedInfo.beam(i).shape{phase}(j).rightLeafPos_I;
                variables.rightLeafPos_F = updatedInfo.beam(i).shape{phase}(j).rightLeafPos_F;
                
                vectorIndices.DAOindex = nnz([updatedInfo.propVMAT.beam(1:i).DAOBeam]);
                vectorIndices.vectorIx_LI = vectorIx_LI;
                vectorIndices.vectorIx_LF = vectorIx_LF;
                vectorIndices.vectorIx_RI = vectorIx_RI;
                vectorIndices.vectorIx_RF = vectorIx_RF;
            else
                
                
                variables.weight_last = updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(j).weight;
                variables.weight_next = updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(j).weight;
                variables.jacobiScale_last = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(1).jacobiScale;
                variables.jacobiScale_next = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(1).jacobiScale;
                variables.leftLeafPos_I = updatedInfo.beam(i).shape{phase}(j).leftLeafPos_I;
                variables.leftLeafPos_F = updatedInfo.beam(i).shape{phase}(j).leftLeafPos_F;
                variables.rightLeafPos_I = updatedInfo.beam(i).shape{phase}(j).rightLeafPos_I;
                variables.rightLeafPos_F = updatedInfo.beam(i).shape{phase}(j).rightLeafPos_F;
                
                variables.time_last = updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).time;
                variables.time_next = updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).time;
                variables.time = updatedInfo.beam(i).time;
                
                variables.fracFromLastOptI = fracFromLastOptI;
                variables.fracFromLastOptF = fracFromLastOptF;
                variables.fracFromNextOptI = fracFromNextOptI;
                variables.fracFromNextOptF = fracFromNextOptF;
                variables.fracFromLastOpt = fracFromLastOpt;
                
                variables.doseAngleBordersDiff = updatedInfo.propVMAT.beam(i).doseAngleBordersDiff;
                variables.doseAngleBordersDiff_last = updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).doseAngleBordersDiff;
                variables.doseAngleBordersDiff_next = updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).doseAngleBordersDiff;
                variables.timeFacCurr_last = updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).timeFacCurr;
                variables.timeFacCurr_next = updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).timeFacCurr;
                variables.fracFromLastDAO = updatedInfo.propVMAT.beam(i).fracFromLastDAO;
                variables.timeFracFromLastDAO = updatedInfo.propVMAT.beam(i).timeFracFromLastDAO;
                variables.timeFracFromNextDAO = updatedInfo.propVMAT.beam(i).timeFracFromNextDAO;
                
                % only part necessary to change for discrete?
                vectorIndices.DAOindex_last = nnz([updatedInfo.propVMAT.beam(1:updatedInfo.propVMAT.beam(i).lastDAOIndex).DAOBeam]);
                vectorIndices.DAOindex_next = nnz([updatedInfo.propVMAT.beam(1:updatedInfo.propVMAT.beam(i).nextDAOIndex).DAOBeam]);
                vectorIndices.vectorIx_LF_last = vectorIx_LF_last;
                vectorIndices.vectorIx_LI_next = vectorIx_LI_next;
                vectorIndices.vectorIx_RF_last = vectorIx_RF_last;
                vectorIndices.vectorIx_RI_next = vectorIx_RI_next;
            end
            
            counters.bixelJApVec_offset = bixelJApVec_offset;
            
            % calculate bixel weight and derivative in function
            [w,bixelJApVec_vec,bixelJApVec_i,bixelJApVec_j,sqrtSumGradSq,shapeMap,counters] = ...
                matRad_bixWeightAndGrad(calcOptions,mlcOptions,variables,vectorIndices,counters,w,bixelJApVec_vec,bixelJApVec_i,bixelJApVec_j);
            
            bixelJApVec_offset = counters.bixelJApVec_offset;
            
            % save the tempMap and sqrtSumGradSq
            updatedInfo.beam(i).shape{phase}(j).shapeMap = shapeMap;
            updatedInfo.beam(i).shape{phase}(j).sqrtSumGradSq = sqrtSumGradSq;
        end
    end
    if updatedInfo.runVMAT
        bixelJApVec_offset = 0;
    end
end

% save bixelWeight, apertureVector, and Jacobian between the two
updatedInfo.bixelWeights = w;
updatedInfo.apertureVector = apertureInfoVect;
updatedInfo.bixelJApVec = cell(apertureInfo.numPhases,1);
for phase = 1:apertureInfo.numPhases
    updatedInfo.bixelJApVec{phase} = sparse(bixelJApVec_i{phase},bixelJApVec_j{phase},bixelJApVec_vec{phase},numel(apertureInfoVect),updatedInfo.totalNumOfBixels);
end

end