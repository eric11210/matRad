function updatedInfo = matRad_daoVec2ApertureInfo_VMATrecalcDynamic(apertureInfo,apertureInfoVect)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to translate the vector representation of the aperture
% shape and weight into an aperture info struct. At the same time, the
% updated bixel weight vector w is computed and a vector listing the
% correspondence between leaf tips and bixel indices for gradient
% calculation
%
% call
%   updatedInfo = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfoVect)
%
% input
%   apertureInfo:     aperture shape info struct
%   apertureInfoVect: aperture weights and shapes parameterized as vector
%   touchingFlag:     if this is one, clean up instances of leaf touching,
%                     otherwise, do not
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
% Copyright 2015 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function to update the apertureInfo struct after the each iteraton of the
% optimization

% initializing variables
updatedInfo = apertureInfo;

updatedInfo.apertureVector = apertureInfoVect;

shapeInd = 1;

% options for bixel and Jacobian calculation
mlcOptions.bixelWidth = apertureInfo.bixelWidth;
calcOptions.continuousAperture = updatedInfo.propVMAT.continuousAperture;
vectorIndices.totalNumOfShapes = apertureInfo.totalNumOfShapes;
vectorIndices.timeOffset = (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases;

w = cell(apertureInfo.numPhases,1);
w(:) = {zeros(apertureInfo.totalNumOfBixels,1)};
% dummy
bixelJApVec_vec = cell(apertureInfo.numPhases,1);
bixelJApVec_i = cell(apertureInfo.numPhases,1);
bixelJApVec_j = cell(apertureInfo.numPhases,1);


%change this to eliminate the first unused entries (which pertain to the
%weights of the aprtures, and to make the bixelIndices work when doing VMAT
%(and we need to potentially interpolate between control points)
indVect = NaN*ones(2*apertureInfo.doseTotalNumOfLeafPairs,1);
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

%Interpolate segment between adjacent optimized gantry angles.
% Include in updatedInfo, but NOT the vector (since these are not
% optimized by DAO).  Also update bixel weights to include these.


%% update the shapeMaps
% here the new colimator positions are used to create new shapeMaps that
% now include decimal values instead of binary

calcOptions.saveJacobian = false;

% loop over all phases
for phase = 1:apertureInfo.numPhases
    % loop over all beams
    for i = 1:numel(updatedInfo.beam)
        
        %posOfRightCornerPixel = apertureInfo.beam(i).posOfCornerBixel(1) + (size(apertureInfo.beam(i).bixelIndMap,2)-1)*apertureInfo.bixelWidth;
        
        % pre compute left and right bixel edges
        edges_l = updatedInfo.beam(i).posOfCornerBixel(1)...
            + ([1:size(apertureInfo.beam(i).bixelIndMap,2)]-1-1/2)*updatedInfo.bixelWidth;
        edges_r = updatedInfo.beam(i).posOfCornerBixel(1)...
            + ([1:size(apertureInfo.beam(i).bixelIndMap,2)]-1+1/2)*updatedInfo.bixelWidth;
        
        % get dimensions of 2d matrices that store shape/bixel information
        n = apertureInfo.beam(i).numOfActiveLeafPairs;
        
        % we are necessarily doing VMAT
        numOfShapes = 1;
        
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
        
        % loop over all shapes
        for j = 1:numOfShapes
            
            % no need to update weights or anything from the vector, just
            % extract the weights and leaf positions from the apertureInfo
            
            weight = updatedInfo.beam(i).shape{phase}(j).weight;
            if isfield(updatedInfo.beam(i).shape{phase}(j),'weight_I')
                weight_I = updatedInfo.beam(i).shape{phase}(j).weight_I;
                weight_F = updatedInfo.beam(i).shape{phase}(j).weight_F;
            else
                %only happens at original angular resolution
                weight_I = weight.*updatedInfo.beam(i).doseAngleBorderCentreDiff(1)./updatedInfo.beam(i).doseAngleBordersDiff;
                weight_F = weight.*updatedInfo.beam(i).doseAngleBorderCentreDiff(2)./updatedInfo.beam(i).doseAngleBordersDiff;
            end
            
            if weight_I+weight_F ~= weight
                %sometimes the sum is different than one by ~10^-16
                %(rounding error in the division)
                weight_F = weight-weight_I;
            end
            
            %% enter in variables and options
            variables.phase = phase;
            counters.bixelJApVec_offset = 0;
            
            %%%%%%%%%%%%%%%%
            %do initial and final arc separately, more accurate
            %calculation
            
            %INITIAL
            variables.weight            = weight_I;
            variables.leftLeafPos_I     = updatedInfo.beam(i).shape{phase}(j).leftLeafPos_I;
            variables.leftLeafPos_F     = updatedInfo.beam(i).shape{phase}(j).leftLeafPos;
            variables.rightLeafPos_I    = updatedInfo.beam(i).shape{phase}(j).rightLeafPos_I;
            variables.rightLeafPos_F    = updatedInfo.beam(i).shape{phase}(j).rightLeafPos;
            
            % calculate bixel weight and derivative in function
            [w,~,bixelJApVec_i,bixelJApVec_j,~,shapeMap_I,counters] = ...
                matRad_bixWeightAndGrad(calcOptions,mlcOptions,variables,vectorIndices,counters,w,bixelJApVec_vec,bixelJApVec_i,bixelJApVec_j);
            
            %FINAL
            variables.weight            = weight_F;
            variables.leftLeafPos_I     = updatedInfo.beam(i).shape{phase}(j).leftLeafPos;
            variables.leftLeafPos_F     = updatedInfo.beam(i).shape{phase}(j).leftLeafPos_F;
            variables.rightLeafPos_I    = updatedInfo.beam(i).shape{phase}(j).rightLeafPos;
            variables.rightLeafPos_F    = updatedInfo.beam(i).shape{phase}(j).rightLeafPos_F;
            
            % calculate bixel weight and derivative in function
            [w,~,bixelJApVec_i,bixelJApVec_j,~,shapeMap_F,counters] = ...
                matRad_bixWeightAndGrad(calcOptions,mlcOptions,variables,vectorIndices,counters,w,bixelJApVec_vec,bixelJApVec_i,bixelJApVec_j);
            
            % save the tempMap
            shapeMap = (shapeMap_I.*weight_I+shapeMap_F.*weight_F)./weight;
            updatedInfo.beam(i).shape{phase}(j).shapeMap = shapeMap;
            
            
            % increment shape index
            shapeInd = shapeInd +1;
        end
    end
    
end

% save bixelWeight, apertureVector, and Jacobian between the two
updatedInfo.bixelWeights = w;
updatedInfo.apertureVector = apertureInfoVect;

end