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
w = cell(apertureInfo.numPhases,1);
w(:) = {zeros(apertureInfo.totalNumOfBixels,1)};

% initializing variables
updatedInfo = apertureInfo;

updatedInfo.apertureVector = apertureInfoVect;

shapeInd = 1;

indVect = NaN*ones(2*apertureInfo.doseTotalNumOfLeafPairs*apertureInfo.numPhases,1);
offset = 0;

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

if updatedInfo.runVMAT && updatedInfo.propVMAT.continuousAperture
    
    % Jacobian matrix to be used in the DAO gradient function
    % this tells us the gradient of a particular bixel with respect to an
    % element in the apertureVector (aperture weight or leaf position)
    % store as a vector for now, convert to sparse matrix later
    bixelJApVec_sz = updatedInfo.totalNumOfOptBixels*5+(updatedInfo.totalNumOfBixels-updatedInfo.totalNumOfOptBixels)*12;
    bixelJApVec_vec = cell(apertureInfo.numPhases,1);
    bixelJApVec_vec(:) = {zeros(1,bixelJApVec_sz)};
    % For optimized beams: 5 = (1 from weights) + (2 from left leaf positions (I and F)) + (2 from
    % right leaf positions (I and F))
    % For interpolated beams: multiply this number times 2 (influenced by the
    % one before and the one after), then add 2 (influenced by the time of the
    % times before and after)
    
    % vector indices
    bixelJApVec_i = cell(apertureInfo.numPhases,1);
    bixelJApVec_i(:) = {zeros(1,bixelJApVec_sz)};
    % bixel indices
    bixelJApVec_j = cell(apertureInfo.numPhases,1);
    bixelJApVec_j(:) = {zeros(1,bixelJApVec_sz)};
    % offset
    bixelJApVec_offset = 0;
end


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
            
            if updatedInfo.propVMAT.continuousAperture
                centres = (edges_l+edges_r)/2;
                widths = edges_r-edges_l;
                numBix = size(apertureInfo.beam(i).bixelIndMap,2);
            end
        else
            numOfShapes = updatedInfo.beam(i).numOfShapes;
        end
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
                    vectorIx_R = vectorIx_L+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
                    leftLeafPos  = apertureInfoVect(vectorIx_L);
                    rightLeafPos = apertureInfoVect(vectorIx_R);
                    
                    % update information in shape structure
                    updatedInfo.beam(i).shape{phase}(j).leftLeafPos  = leftLeafPos;
                    updatedInfo.beam(i).shape{phase}(j).rightLeafPos = rightLeafPos;
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
                    % obtain leaf positions at last DAO beam
                    vectorIx_L_last = updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(j).vectorOffset + ([1:n]-1);
                    vectorIx_R_last = vectorIx_L_last+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
                    leftLeafPos_last = apertureInfoVect(vectorIx_L_last);
                    rightLeafPos_last = apertureInfoVect(vectorIx_R_last);
                    
                    % obtain leaf positions at next DAO beam
                    vectorIx_L_next = updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(j).vectorOffset + ([1:n]-1);
                    vectorIx_R_next = vectorIx_L_next+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
                    leftLeafPos_next = apertureInfoVect(vectorIx_L_next);
                    rightLeafPos_next = apertureInfoVect(vectorIx_R_next);
                    
                    % interpolate leaf positions
                    leftLeafPos = updatedInfo.propVMAT.beam(i).fracFromLastDAO*leftLeafPos_last+(1-updatedInfo.propVMAT.beam(i).fracFromLastDAO)*leftLeafPos_next;
                    rightLeafPos = updatedInfo.propVMAT.beam(i).fracFromLastDAO*rightLeafPos_last+(1-updatedInfo.propVMAT.beam(i).fracFromLastDAO)*rightLeafPos_next;
                    
                    % update information in shape structure
                    updatedInfo.beam(i).shape{phase}(j).leftLeafPos  = leftLeafPos;
                    updatedInfo.beam(i).shape{phase}(j).rightLeafPos = rightLeafPos;
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
                w{phase}(currBixelIx) = w{phase}(currBixelIx) + tempMap(tempMapIx)*updatedInfo.beam(i).shape{phase}(j).weight;
                
                % save the tempMap (we need to apply a positivity operator !)
                updatedInfo.beam(i).shape{phase}(j).shapeMap = (tempMap  + abs(tempMap))  / 2;
                
                if ~updatedInfo.runVMAT || updatedInfo.propVMAT.beam(i).DAOBeam
                    % increment shape index
                    shapeInd = shapeInd +1;
                end
            else
                %% continuous aperture
                
                % calculate bixel weight and derivative in function
                
                % the weight here is going to depend on the beginning and
                % ending phases
                [w,bixelJApVec_vec,bixelJApVec_i,bixelJApVec_j,tempMap] = ...
                    matRad_bixWeightAndGrad('continuous',updatedInfo.propVMAT.beam(i).DAOBeam, ...
                    apertureInfo.beam(i).lim_l,apertureInfo.beam(i).lim_r,edges_l,edges_r,centres,widths,n,numBix,apertureInfo.bixelWidth, ...
                    apertureInfo.beam(i).bixelIndMap,updatedInfo.beam(i).shape{phase}(j).weight,apertureInfo.beam(i).shape{phase}(1).jacobiScale, ...
                    w,bixelJApVec_vec,bixelJApVec_i,bixelJApVec_j, ...
                    updatedInfo.beam(i).shape{phase}(j).leftLeafPos_I,updatedInfo.beam(i).shape{phase}(j).leftLeafPos_F, ...
                    updatedInfo.beam(i).shape{phase}(j).rightLeafPos_I,updatedInfo.beam(i).shape{phase}(j).rightLeafPos_F, ...
                    vectorIx_LI,vectorIx_LF,vectorIx_RI,vectorIx_RF);
                
                % save the tempMap
                updatedInfo.beam(i).shape{phase}(j).shapeMap = tempMap;
                
                if updatedInfo.propVMAT.beam(i).DAOBeam
                    % indices
                    numSaveBixel = nnz(saveBixelIx);
                    vectorIx_LI = repmat(vectorIx_LI',1,numBix);
                    vectorIx_LF = repmat(vectorIx_LF',1,numBix);
                    vectorIx_RI = repmat(vectorIx_RI',1,numBix);
                    vectorIx_RF = repmat(vectorIx_RF',1,numBix);
                    
                    % wrt weight
                    bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = updatedInfo.beam(i).shape{phase}(j).shapeMap(saveBixelIx)./apertureInfo.beam(i).shape{phase}(1).jacobiScale;
                    bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = shapeInd;
                    bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
                    bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
                    
                    % wrt initial left
                    bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = dUl_dLI(saveBixelIx)*updatedInfo.beam(i).shape{phase}(j).weight;
                    bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_LI(saveBixelIx);
                    bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
                    bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
                    
                    % wrt final left
                    bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = dUl_dLF(saveBixelIx)*updatedInfo.beam(i).shape{phase}(j).weight;
                    bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_LF(saveBixelIx);
                    bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
                    bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
                    
                    % wrt initial right
                    bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = -dCr_dRI(saveBixelIx)*updatedInfo.beam(i).shape{phase}(j).weight;
                    bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_RI(saveBixelIx);
                    bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
                    bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
                    
                    % wrt final right
                    bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = -dCr_dRF(saveBixelIx)*updatedInfo.beam(i).shape{phase}(j).weight;
                    bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_RF(saveBixelIx);
                    bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
                    bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
                    
                    % store information for Jacobi preconditioning
                    updatedInfo.beam(i).shape{phase}(j).sqrtSumGradSq = sqrt(mean([sum(dUl_dLI.^2,2); sum(dUl_dLF.^2,2); sum(dCr_dRI.^2,2); sum(dCr_dRF.^2,2)]));
                    
                    % increment shape index
                    shapeInd = shapeInd +1;
                else
                    % indices
                    saveBixelIx = ~isnan(apertureInfo.beam(i).bixelIndMap);
                    numSaveBixel = nnz(saveBixelIx);
                    vectorIx_LF_last = repmat(vectorIx_LF_last',1,numBix);
                    vectorIx_LI_next = repmat(vectorIx_LI_next',1,numBix);
                    vectorIx_RF_last = repmat(vectorIx_RF_last',1,numBix);
                    vectorIx_RI_next = repmat(vectorIx_RI_next',1,numBix);
                    
                    % leaf interpolation fractions/weights
                    fracFromLastOptI = repmat(fracFromLastOptI,1,numBix);
                    fracFromLastOptF = repmat(fracFromLastOptF,1,numBix);
                    fracFromNextOptI = repmat(fracFromNextOptI,1,numBix);
                    fracFromNextOptF = repmat(fracFromNextOptF,1,numBix);
                    
                    % wrt last weight
                    bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = fracFromLastOpt*(updatedInfo.beam(i).time./updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).time)*updatedInfo.beam(i).shape{phase}(j).shapeMap(saveBixelIx) ...
                        ./apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(1).jacobiScale;
                    %bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = (updatedInfo.beam(i).doseAngleBordersDiff*fracFromLastOpt*updatedInfo.beam(apertureInfo.beam(i).lastOptIndex).gantryRot ...
                    %/(updatedInfo.beam(apertureInfo.beam(i).lastOptIndex).doseAngleBordersDiff*updatedInfo.beam(i).gantryRot))*updatedInfo.beam(i).shape(j).shapeMap(saveBixelIx) ...
                    %./ apertureInfo.beam(apertureInfo.beam(i).lastOptIndex).shape(1).jacobiScale;
                    bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = nnz([updatedInfo.propVMAT.beam(1:updatedInfo.propVMAT.beam(i).lastDAOIndex).DAOBeam])+(phase-1)*apertureInfo.totalNumOfShapes;
                    bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
                    bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
                    
                    % wrt next weight
                    bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = (1-fracFromLastOpt)*(updatedInfo.beam(i).time./updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).time)*updatedInfo.beam(i).shape{phase}(j).shapeMap(saveBixelIx) ...
                        ./apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(1).jacobiScale;
                    %bixelJApVec_vec(bixelJApVec_offset+(1:numSaveBixel)) = (updatedInfo.beam(i).doseAngleBordersDiff*(1-fracFromLastOpt)*updatedInfo.beam(apertureInfo.beam(i).nextOptIndex).gantryRot ...
                    %/(updatedInfo.beam(apertureInfo.beam(i).nextOptIndex).doseAngleBordersDiff*updatedInfo.beam(i).gantryRot))*updatedInfo.beam(i).shape(j).shapeMap(saveBixelIx) ...
                    %./ apertureInfo.beam(apertureInfo.beam(i).nextOptIndex).shape(1).jacobiScale;
                    bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = nnz([updatedInfo.propVMAT.beam(1:updatedInfo.propVMAT.beam(i).nextDAOIndex).DAOBeam])+(phase-1)*apertureInfo.totalNumOfShapes;
                    bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
                    bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
                    
                    
                    %updatedInfo.beam(i).shape(j).leftLeafPos_I = fracFromLastOptI.*leftLeafPos_F_last+fracFromNextOptI.*leftLeafPos_I_next;
                    %updatedInfo.beam(i).shape(j).rightLeafPos_I = fracFromLastOptI.*rightLeafPos_F_last+fracFromNextOptI.*rightLeafPos_I_next;
                    
                    %updatedInfo.beam(i).shape(j).leftLeafPos_F = fracFromLastOptF.*leftLeafPos_F_last+fracFromNextOptF.*leftLeafPos_I_next;
                    %updatedInfo.beam(i).shape(j).rightLeafPos_F = fracFromLastOptF.*rightLeafPos_F_last+fracFromNextOptF.*rightLeafPos_I_next;
                    
                    % wrt initial left (optimization vector)
                    % initial (interpolated arc)
                    bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = fracFromLastOptI(saveBixelIx).*dUl_dLI(saveBixelIx)*updatedInfo.beam(i).shape{phase}(j).weight;
                    bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_LF_last(saveBixelIx);
                    bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
                    bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
                    % final (interpolated arc)
                    bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = fracFromLastOptF(saveBixelIx).*dUl_dLF(saveBixelIx)*updatedInfo.beam(i).shape{phase}(j).weight;
                    bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_LF_last(saveBixelIx);
                    bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
                    bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
                    
                    % wrt final left (optimization vector)
                    % initial (interpolated arc)
                    bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = fracFromNextOptI(saveBixelIx).*dUl_dLI(saveBixelIx)*updatedInfo.beam(i).shape{phase}(j).weight;
                    bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_LI_next(saveBixelIx);
                    bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
                    bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
                    % final (interpolated arc)
                    bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = fracFromNextOptF(saveBixelIx).*dUl_dLF(saveBixelIx)*updatedInfo.beam(i).shape{phase}(j).weight;
                    bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_LI_next(saveBixelIx);
                    bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
                    bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
                    
                    % wrt initial right (optimization vector)
                    % initial (interpolated arc)
                    bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = -fracFromLastOptI(saveBixelIx).*dCr_dRI(saveBixelIx)*updatedInfo.beam(i).shape{phase}(j).weight;
                    bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_RF_last(saveBixelIx);
                    bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
                    bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
                    % final (interpolated arc)
                    bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = -fracFromLastOptF(saveBixelIx).*dCr_dRF(saveBixelIx)*updatedInfo.beam(i).shape{phase}(j).weight;
                    bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_RF_last(saveBixelIx);
                    bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
                    bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
                    
                    % wrt final right (optimization vector)
                    % initial (interpolated arc)
                    bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = -fracFromNextOptI(saveBixelIx).*dCr_dRI(saveBixelIx)*updatedInfo.beam(i).shape{phase}(j).weight;
                    bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_RI_next(saveBixelIx);
                    bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
                    bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
                    % final (interpolated arc)
                    bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = -fracFromNextOptF(saveBixelIx).*dCr_dRF(saveBixelIx)*updatedInfo.beam(i).shape{phase}(j).weight;
                    bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = vectorIx_RI_next(saveBixelIx);
                    bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
                    bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
                    
                    % wrt last time
                    bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = (updatedInfo.propVMAT.beam(i).doseAngleBordersDiff.*updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).timeFacCurr) ...
                        .*(-updatedInfo.propVMAT.beam(i).fracFromLastDAO.*updatedInfo.propVMAT.beam(i).timeFracFromNextDAO.*(updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(j).weight./updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).doseAngleBordersDiff).*(updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).time./updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).time.^2) ...
                        +(1-updatedInfo.propVMAT.beam(i).fracFromLastDAO).*updatedInfo.propVMAT.beam(i).timeFracFromLastDAO.*(updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(j).weight./updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).doseAngleBordersDiff).*(1./updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).time)) ...
                        * updatedInfo.beam(i).shape{phase}(j).shapeMap(saveBixelIx);
                    bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = nnz([updatedInfo.propVMAT.beam(1:updatedInfo.propVMAT.beam(i).lastDAOIndex).DAOBeam])+(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases;
                    bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
                    bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
                    
                    % wrt next time
                    bixelJApVec_vec{phase}(bixelJApVec_offset+(1:numSaveBixel)) = (updatedInfo.propVMAT.beam(i).doseAngleBordersDiff.*updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).timeFacCurr) ...
                        .*(updatedInfo.propVMAT.beam(i).fracFromLastDAO.*updatedInfo.propVMAT.beam(i).timeFracFromNextDAO.*(updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(j).weight./updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).doseAngleBordersDiff).*(1./updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).time) ...
                        -(1-updatedInfo.propVMAT.beam(i).fracFromLastDAO).*updatedInfo.propVMAT.beam(i).timeFracFromLastDAO.*(updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(j).weight./updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).doseAngleBordersDiff).*(updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).time./updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).time.^2)) ...
                        * updatedInfo.beam(i).shape{phase}(j).shapeMap(saveBixelIx);
                    bixelJApVec_i{phase}(bixelJApVec_offset+(1:numSaveBixel)) = nnz([updatedInfo.propVMAT.beam(1:updatedInfo.propVMAT.beam(i).nextDAOIndex).DAOBeam])+(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases;
                    bixelJApVec_j{phase}(bixelJApVec_offset+(1:numSaveBixel)) = apertureInfo.beam(i).bixelIndMap(saveBixelIx);
                    bixelJApVec_offset = bixelJApVec_offset+numSaveBixel;
                end
            end
        end
    end
    if updatedInfo.runVMAT && updatedInfo.propVMAT.continuousAperture
        bixelJApVec_offset = 0;
    end
end


updatedInfo.bixelWeights = w;
updatedInfo.apertureVector = apertureInfoVect;
if ~updatedInfo.propVMAT.continuousAperture
    updatedInfo.bixelIndices = indVect;
else
    updatedInfo.bixelJApVec = cell(apertureInfo.numPhases,1);
    for phase = 1:apertureInfo.numPhases
        updatedInfo.bixelJApVec{phase} = sparse(bixelJApVec_i{phase},bixelJApVec_j{phase},bixelJApVec_vec{phase},numel(apertureInfoVect),updatedInfo.totalNumOfBixels,bixelJApVec_sz);
    end
end

end