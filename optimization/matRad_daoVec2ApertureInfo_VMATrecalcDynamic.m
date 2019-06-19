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

% options for bixel and Jacobian calculation
static.bixelWidth = apertureInfo.bixelWidth;
variable.totalNumOfShapes = apertureInfo.totalNumOfShapes;

w = cell(apertureInfo.numPhases,1);
w(:) = {zeros(apertureInfo.totalNumOfBixels,1)};
%{
if updatedInfo.runVMAT && ~all([updatedInfo.propVMAT.beam.DAOBeam])
    j = 1;
    for phase = 1:apertureInfo.numPhases
        for i = 1:numel(updatedInfo.beam)
            if updatedInfo.propVMAT.beam(i).DAOBeam
                % update the shape weight
                % rescale the weight from the vector using the previous
                % iteration scaling factor
                updatedInfo.beam(i).shape{phase}(j).weight = apertureInfoVect(updatedInfo.beam(i).shape{phase}(j).weightOffset)./updatedInfo.beam(i).shape{phase}(j).jacobiScale;
                
                updatedInfo.beam(i).shape{phase}(j).MU = updatedInfo.beam(i).shape{phase}(j).weight*updatedInfo.weightToMU;
                if phase == 1
                    updatedInfo.beam(i).time = apertureInfoVect((updatedInfo.totalNumOfShapes+updatedInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases+updatedInfo.propVMAT.beam(i).DAOIndex)*updatedInfo.propVMAT.beam(i).timeFacCurr;
                    updatedInfo.beam(i).gantryRot = updatedInfo.propVMAT.beam(i).doseAngleBordersDiff/updatedInfo.beam(i).time;
                end
                updatedInfo.beam(i).shape{phase}(j).MURate = updatedInfo.beam(i).shape{phase}(j).MU./updatedInfo.beam(i).time;
            end
        end
    end
end
%}

%% update the shapeMaps
% here the new colimator positions are used to create new shapeMaps that
% now include decimal values instead of binary

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
    
    %weightFactor_I = updatedInfo.propVMAT.beam(i).doseAngleBorderCentreDiff(1)./updatedInfo.propVMAT.beam(i).doseAngleBordersDiff;
    %weightFactor_F = updatedInfo.propVMAT.beam(i).doseAngleBorderCentreDiff(2)./updatedInfo.propVMAT.beam(i).doseAngleBordersDiff;
    
    % we are necessarily doing VMAT
    numOfShapes = 1;
    
    static.lim_r        = apertureInfo.beam(i).lim_r;
    static.edges_l      = edges_l;
    static.edges_r      = edges_r;
    static.centres      = (edges_l+edges_r)/2;
    static.widths       = edges_r-edges_l;
    static.numRow       = n;
    static.numCol       = size(apertureInfo.beam(i).bixelIndMap,2);
    static.numBix       = apertureInfo.beam(i).numBixels;
    static.bixelIndMap  = apertureInfo.beam(i).bixelIndMap-max(apertureInfo.beam(i).bixelIndMap(:))+apertureInfo.beam(i).numBixels;
    static.DAOBeam      = updatedInfo.propVMAT.beam(i).DAOBeam;
    static.bixIndVec    = 1:static.numBix;
    
    % dummies
    variable.probability_dTVec          = ones(numel(apertureInfo.beam),1);
    variable.tIx_Vec                    = ones(apertureInfo.totalNumOfShapes,1);
    variable.numShapbixIndVec           = 1:(variable.totalNumOfShapes*static.numBix);
    variable.DAOindex                   = 1;
    variable.jacobiScale                = 1;
    variable.vectorIx_LI                = ones(static.numRow,1);
    variable.vectorIx_LF                = ones(static.numRow,1);
    variable.vectorIx_RI                = ones(static.numRow,1);
    variable.vectorIx_RF                = ones(static.numRow,1);
    variable.bixelJApVec_sz             = static.numBix.*(7+variable.totalNumOfShapes);
    variable.arcF                       = false;
    
    % loop over all shapes
    for j = 1:numOfShapes
        
        % loop over all phases
        for phase = 1:apertureInfo.numPhases
            
            
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
            
            %%%%%%%%%%%%%%%%
            %do initial and final arc separately, more accurate
            %calculation
            
            %INITIAL
            variable.weight             = weight_I;
            variable.weightFactor_I     = 1;
            variable.weightFactor_F     = 1;
            variable.probability        = apertureInfo.motionModel.probPhase(phase);
            
            if updatedInfo.propVMAT.continuousAperture
                leftLeafPosI    = updatedInfo.beam(i).shape{phase}(j).leftLeafPos_I;
                leftLeafPosF    = updatedInfo.beam(i).shape{phase}(j).leftLeafPos;
                rightLeafPosI   = updatedInfo.beam(i).shape{phase}(j).rightLeafPos_I;
                rightLeafPosF   = updatedInfo.beam(i).shape{phase}(j).rightLeafPos;
            else
                leftLeafPosI    = updatedInfo.beam(i).shape{phase}(j).leftLeafPos;
                leftLeafPosF    = updatedInfo.beam(i).shape{phase}(j).leftLeafPos;
                rightLeafPosI   = updatedInfo.beam(i).shape{phase}(j).rightLeafPos;
                rightLeafPosF   = updatedInfo.beam(i).shape{phase}(j).rightLeafPos;
            end
            
            % set the initial leaf positions to the minimum leaf positions
            % always, instead of the leaf positions at the actual beginning
            % of the arc
            % this simplifies the calculation
            variable.leftLeafPosI   = min([leftLeafPosI,leftLeafPosF],[],2);
            variable.leftLeafPosF   = max([leftLeafPosI,leftLeafPosF],[],2);
            variable.rightLeafPosI  = min([rightLeafPosI,rightLeafPosF],[],2);
            variable.rightLeafPosF  = max([rightLeafPosI,rightLeafPosF],[],2);
            
            % find bixel indices where leaves are located
            variable.xPosIndLeftLeafI    = min(floor((variable.leftLeafPosI-edges_l(1))./apertureInfo.bixelWidth)+1,static.numCol);
            variable.xPosIndLeftLeafF    = max(ceil((variable.leftLeafPosF-edges_r(1))./apertureInfo.bixelWidth)+1,1);
            variable.xPosIndRightLeafI   = min(floor((variable.rightLeafPosI-edges_l(1))./apertureInfo.bixelWidth)+1,static.numCol);
            variable.xPosIndRightLeafF   = max(ceil((variable.rightLeafPosF-edges_r(1))./apertureInfo.bixelWidth)+1,1);
            %
            variable.xPosLinearIndLeftLeafI      = sub2ind([static.numRow static.numCol],(1:static.numRow)',variable.xPosIndLeftLeafI);
            variable.xPosLinearIndLeftLeafF      = sub2ind([static.numRow static.numCol],(1:static.numRow)',variable.xPosIndLeftLeafF);
            variable.xPosLinearIndRightLeafI     = sub2ind([static.numRow static.numCol],(1:static.numRow)',variable.xPosIndRightLeafI);
            variable.xPosLinearIndRightLeafF     = sub2ind([static.numRow static.numCol],(1:static.numRow)',variable.xPosIndRightLeafF);
            
            % calculate bixel weight and derivative in function
            tempResults = matRad_bixWeightAndGrad(static,variable);
            
            % put tempResults into results
            w{phase}(apertureInfo.beam(i).bixelIndMap(~isnan(apertureInfo.beam(i).bixelIndMap))) = w{phase}(apertureInfo.beam(i).bixelIndMap(~isnan(apertureInfo.beam(i).bixelIndMap)))+tempResults.w;
            shapeMap_I = tempResults.shapeMap;
            
            
            %FINAL
            variable.weight             = weight_F;
            variable.weightFactor_I     = 1;
            variable.weightFactor_F     = 1;
            variable.probability        = apertureInfo.motionModel.probPhase(phase);
            if updatedInfo.propVMAT.continuousAperture
                leftLeafPosI    = updatedInfo.beam(i).shape{phase}(j).leftLeafPos;
                leftLeafPosF    = updatedInfo.beam(i).shape{phase}(j).leftLeafPos_F;
                rightLeafPosI   = updatedInfo.beam(i).shape{phase}(j).rightLeafPos;
                rightLeafPosF   = updatedInfo.beam(i).shape{phase}(j).rightLeafPos_F;
            else
                leftLeafPosI    = updatedInfo.beam(i).shape{phase}(j).leftLeafPos;
                leftLeafPosF    = updatedInfo.beam(i).shape{phase}(j).leftLeafPos;
                rightLeafPosI   = updatedInfo.beam(i).shape{phase}(j).rightLeafPos;
                rightLeafPosF   = updatedInfo.beam(i).shape{phase}(j).rightLeafPos;
            end
            
            % set the initial leaf positions to the minimum leaf positions
            % always, instead of the leaf positions at the actual beginning
            % of the arc
            % this simplifies the calculation
            variable.leftLeafPosI   = min([leftLeafPosI,leftLeafPosF],[],2);
            variable.leftLeafPosF   = max([leftLeafPosI,leftLeafPosF],[],2);
            variable.rightLeafPosI  = min([rightLeafPosI,rightLeafPosF],[],2);
            variable.rightLeafPosF  = max([rightLeafPosI,rightLeafPosF],[],2);
            
            % find bixel indices where leaves are located
            variable.xPosIndLeftLeafI    = min(floor((variable.leftLeafPosI-edges_l(1))./apertureInfo.bixelWidth)+1,static.numCol);
            variable.xPosIndLeftLeafF    = max(ceil((variable.leftLeafPosF-edges_r(1))./apertureInfo.bixelWidth)+1,1);
            variable.xPosIndRightLeafI   = min(floor((variable.rightLeafPosI-edges_l(1))./apertureInfo.bixelWidth)+1,static.numCol);
            variable.xPosIndRightLeafF   = max(ceil((variable.rightLeafPosF-edges_r(1))./apertureInfo.bixelWidth)+1,1);
            %
            variable.xPosLinearIndLeftLeafI      = sub2ind([static.numRow static.numCol],(1:static.numRow)',variable.xPosIndLeftLeafI);
            variable.xPosLinearIndLeftLeafF      = sub2ind([static.numRow static.numCol],(1:static.numRow)',variable.xPosIndLeftLeafF);
            variable.xPosLinearIndRightLeafI     = sub2ind([static.numRow static.numCol],(1:static.numRow)',variable.xPosIndRightLeafI);
            variable.xPosLinearIndRightLeafF     = sub2ind([static.numRow static.numCol],(1:static.numRow)',variable.xPosIndRightLeafF);
            
            % calculate bixel weight and derivative in function
            tempResults = matRad_bixWeightAndGrad(static,variable);
            
            % put tempResults into results
            w{phase}(apertureInfo.beam(i).bixelIndMap(~isnan(apertureInfo.beam(i).bixelIndMap))) = w{phase}(apertureInfo.beam(i).bixelIndMap(~isnan(apertureInfo.beam(i).bixelIndMap)))+tempResults.w;
            shapeMap_F = tempResults.shapeMap;
            
            % save the tempMap
            shapeMap = shapeMap_I+shapeMap_F;
            updatedInfo.beam(i).shape{phase}(j).shapeMap = shapeMap;
            
        end
    end
end

% save bixelWeight, apertureVector, and Jacobian between the two
updatedInfo.bixelWeights = w;
updatedInfo.apertureVector = apertureInfoVect;

end