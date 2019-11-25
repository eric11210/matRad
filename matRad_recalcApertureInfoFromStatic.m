function recalc = matRad_recalcApertureInfoFromStatic(recalc,apertureInfoOld)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to apertures for a different dose resolution
%
% call
%   recalc = matRad_recalcApertureInfo(recalc,apertureInfo)
%
% input
%   recalc:
%   apertureInfo:
%
% output
%   recalc:
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

pln = recalc.pln;
stf = recalc.stf;

apertureInfoNew = apertureInfoOld;
apertureInfoNew = rmfield(apertureInfoNew,'beam');

apertureInfoNew.totalNumOfBixels = sum([stf(:).totalNumOfBixels]);

shapeInd = 1;

if recalc.interpNew
    oldGantryAngles = zeros(1,numel(apertureInfoOld.beam));
    oldLeftLeafPoss = zeros(apertureInfoOld.beam(1).numOfActiveLeafPairs,numel(apertureInfoOld.beam));
    oldRightLeafPoss = zeros(apertureInfoOld.beam(1).numOfActiveLeafPairs,numel(apertureInfoOld.beam));
    for phaseOld = 1:apertureInfoOld.numPhases
        for iOld = 1:numel(apertureInfoOld.beam)
            oldGantryAngles(iOld) = apertureInfoOld.beam(iOld).gantryAngle;
            oldLeftLeafPoss(:,iOld) = apertureInfoOld.beam(iOld).shape{phaseOld}(1).leftLeafPos;
            oldRightLeafPoss(:,iOld) = apertureInfoOld.beam(iOld).shape{phaseOld}(1).rightLeafPos;
        end
    end
end

% MLC parameters:
bixelWidth = stf(1).bixelWidth; % [mm]
numOfMLCLeafPairs = 80;
%     define central leaf pair (here we want the 0mm position to be in the
%     center of a leaf pair (e.g. leaf 41 stretches from -2.5mm to 2.5mm
%     for a bixel/leafWidth of 5mm and 81 leaf pairs)
centralLeafPair = ceil(numOfMLCLeafPairs/2);

% initializing variables
totalNumOfShapes = numel(stf);
% loop over all phases
for phaseOld = 1:apertureInfoOld.numPhases
    for iOld = 1:numel(apertureInfoOld.beam)
        newInd = (apertureInfoOld.propVMAT.beam(iOld).doseAngleBorders(1) <= [stf.gantryAngle] & [stf.gantryAngle] <= apertureInfoOld.propVMAT.beam(iOld).doseAngleBorders(2)).*(1:numel([stf.gantryAngle]));
        newInd(newInd == 0) = [];
        
        totalAmountOfOldWeight = 0;
        
        for iNew = newInd
            % get x- and z-coordinates of bixels
            rayPos_bev = reshape([stf(iNew).ray.rayPos_bev],3,[]);
            X = rayPos_bev(1,:)';
            Z = rayPos_bev(3,:)';
            
            % create ray-map
            maxX = max(X); minX = min(X);
            maxZ = max(Z); minZ = min(Z);
            
            dimX = (maxX-minX)/stf(iNew).bixelWidth + 1;
            dimZ = (maxZ-minZ)/stf(iNew).bixelWidth + 1;
            
            rayMap = zeros(dimZ,dimX);
            
            % get indices for x and z positions
            xPos = (X-minX)/stf(iNew).bixelWidth + 1;
            zPos = (Z-minZ)/stf(iNew).bixelWidth + 1;
            
            % get indices in the ray-map
            indInRay = zPos + (xPos-1)*dimZ;
            
            % fill ray-map
            rayMap(indInRay) = 1;
            
            % create map of bixel indices
            bixelIndMap = NaN * ones(dimZ,dimX);
            bixelIndMap(indInRay) = [1:stf(iNew).numOfRays] + (iNew-1)*stf(1).numOfRays;
            
            % store physical position of first entry in bixelIndMap
            posOfCornerBixel = [minX 0 minZ];
            
            % get leaf limits from the leaf map
            lim_l = NaN * ones(dimZ,1);
            lim_r = NaN * ones(dimZ,1);
            % looping oder leaf pairs
            for l = 1:dimZ
                lim_lInd = find(rayMap(l,:),1,'first');
                lim_rInd = find(rayMap(l,:),1,'last');
                % the physical position [mm] can be calculated from the indices
                lim_l(l) = (lim_lInd-1)*bixelWidth + minX - 1/2*bixelWidth;
                lim_r(l) = (lim_rInd-1)*bixelWidth + minX + 1/2*bixelWidth;
            end
            
            leafPairPos = unique(Z);
            
            % find upmost and downmost leaf pair
            topLeafPairPos = maxZ;
            bottomLeafPairPos = minZ;
            
            topLeafPair = centralLeafPair - topLeafPairPos/bixelWidth;
            bottomLeafPair = centralLeafPair - bottomLeafPairPos/bixelWidth;
            
            % create bool map of active leaf pairs
            isActiveLeafPair = zeros(numOfMLCLeafPairs,1);
            isActiveLeafPair(topLeafPair:bottomLeafPair) = 1;
            
            MLCWindow = [minX-bixelWidth/2 maxX+bixelWidth/2 ...
                minZ-bixelWidth/2 maxZ+bixelWidth/2];
            
            
            % save data for each beam
            apertureInfoNew.beam(iNew).numOfActiveLeafPairs = dimZ;
            apertureInfoNew.beam(iNew).leafPairPos = leafPairPos;
            apertureInfoNew.beam(iNew).isActiveLeafPair = isActiveLeafPair;
            apertureInfoNew.beam(iNew).centralLeafPair = centralLeafPair;
            apertureInfoNew.beam(iNew).lim_l = lim_l;
            apertureInfoNew.beam(iNew).lim_r = lim_r;
            apertureInfoNew.beam(iNew).bixelIndMap = bixelIndMap;
            apertureInfoNew.beam(iNew).posOfCornerBixel = posOfCornerBixel;
            apertureInfoNew.beam(iNew).MLCWindow = MLCWindow;
            apertureInfoNew.beam(iNew).bixOffset = 1+(iNew-1)*dimZ;
            apertureInfoNew.beam(iNew).numBixels = nnz(~isnan(bixelIndMap));
            apertureInfoNew.beam(iNew).shape{phaseOld}(1).vectorOffset = totalNumOfShapes+1+(iNew-1)*dimZ;
            
            %inherit from old beam
            apertureInfoNew.propVMAT.beam(iNew).leafDir = apertureInfoOld.propVMAT.beam(iOld).leafDir;
            
            %specific to new beam
            apertureInfoNew.beam(iNew).gantryAngle = pln.propStf.gantryAngles(iNew);
            apertureInfoNew.propVMAT.beam(iNew).doseAngleBorders = stf(iNew).propVMAT.doseAngleBorders;
            apertureInfoNew.propVMAT.beam(iNew).doseAngleBorderCentreDiff = stf(iNew).propVMAT.doseAngleBorderCentreDiff;
            apertureInfoNew.propVMAT.beam(iNew).doseAngleBordersDiff = stf(iNew).propVMAT.doseAngleBordersDiff;
            apertureInfoNew.propVMAT.beam(iNew).lastDAOIndex = stf(iNew).propVMAT.lastDAOIndex;
            apertureInfoNew.propVMAT.beam(iNew).nextDAOIndex = stf(iNew).propVMAT.lastDAOIndex;
            
            
            amountOfOldSpeed = (min(apertureInfoNew.propVMAT.beam(iNew).doseAngleBorders(2),apertureInfoOld.propVMAT.beam(iOld).doseAngleBorders(2))-max(apertureInfoNew.propVMAT.beam(iNew).doseAngleBorders(1),apertureInfoOld.propVMAT.beam(iOld).doseAngleBorders(1)))./apertureInfoNew.propVMAT.beam(iNew).doseAngleBordersDiff;
            amountOfOldWeight = (min(apertureInfoNew.propVMAT.beam(iNew).doseAngleBorders(2),apertureInfoOld.propVMAT.beam(iOld).doseAngleBorders(2))-max(apertureInfoNew.propVMAT.beam(iNew).doseAngleBorders(1),apertureInfoOld.propVMAT.beam(iOld).doseAngleBorders(1)))./apertureInfoOld.propVMAT.beam(iOld).doseAngleBordersDiff;
            
            totalAmountOfOldWeight = totalAmountOfOldWeight+amountOfOldWeight;
            
            amountOfOldWeight_I = (min(apertureInfoNew.beam(iNew).gantryAngle,apertureInfoOld.propVMAT.beam(iOld).doseAngleBorders(2))-max(apertureInfoNew.propVMAT.beam(iNew).doseAngleBorders(1),apertureInfoOld.propVMAT.beam(iOld).doseAngleBorders(1)))./apertureInfoOld.propVMAT.beam(iOld).doseAngleBordersDiff;
            amountOfOldWeight_F = (min(apertureInfoNew.propVMAT.beam(iNew).doseAngleBorders(2),apertureInfoOld.propVMAT.beam(iOld).doseAngleBorders(2))-max(apertureInfoNew.beam(iNew).gantryAngle,apertureInfoOld.propVMAT.beam(iOld).doseAngleBorders(1)))./apertureInfoOld.propVMAT.beam(iOld).doseAngleBordersDiff;
            
            if ~isfield(apertureInfoNew.beam(iNew),'gantryRot') || isempty(apertureInfoNew.beam(iNew).gantryRot)
                apertureInfoNew.beam(iNew).gantryRot = 0;
                apertureInfoNew.beam(iNew).shape{phaseOld}(1).weight = 0;
                apertureInfoNew.beam(iNew).shape{phaseOld}(1).weight_I = 0;
                apertureInfoNew.beam(iNew).shape{phaseOld}(1).weight_F = 0;
            end
            apertureInfoNew.beam(iNew).gantryRot = amountOfOldSpeed*apertureInfoOld.beam(iOld).gantryRot+apertureInfoNew.beam(iNew).gantryRot;
            
            %recalculate weight, MU
            apertureInfoNew.beam(iNew).shape{phaseOld}(1).weight = apertureInfoNew.beam(iNew).shape{phaseOld}(1).weight+amountOfOldWeight*apertureInfoOld.beam(iOld).shape{phaseOld}(1).weight;
            apertureInfoNew.beam(iNew).shape{phaseOld}(1).weight_I = apertureInfoNew.beam(iNew).shape{phaseOld}(1).weight_I+amountOfOldWeight_I*apertureInfoOld.beam(iOld).shape{phaseOld}(1).weight;
            apertureInfoNew.beam(iNew).shape{phaseOld}(1).weight_F = apertureInfoNew.beam(iNew).shape{phaseOld}(1).weight_F+amountOfOldWeight_F*apertureInfoOld.beam(iOld).shape{phaseOld}(1).weight;
            apertureInfoNew.beam(iNew).MU = apertureInfoNew.beam(iNew).shape{phaseOld}(1).weight.*apertureInfoNew.weightToMU;
            
            apertureInfoNew.beam(iNew).MURate = apertureInfoNew.beam(iNew).MU.*apertureInfoNew.beam(iNew).gantryRot./apertureInfoNew.propVMAT.beam(iNew).doseAngleBordersDiff;
            
            %apertureInfoNew.beam(j).shape{phase}(1).jacobiScale = apertureInfoOld.beam(i).shape{phase}(1).jacobiScale;
            apertureInfoNew.jacobiScale(iNew) = 1;
            apertureInfoNew.beam(iNew).shape{phaseOld}(1).jacobiScale = 1;
            
            if recalc.interpNew
                %interpolate new apertures now so that weights are not
                %overwritten
                apertureInfoNew.beam(iNew).shape{phaseOld}(1).leftLeafPos = (interp1(oldGantryAngles',oldLeftLeafPoss',apertureInfoNew.beam(iNew).gantryAngle))';
                apertureInfoNew.beam(iNew).shape{phaseOld}(1).rightLeafPos = (interp1(oldGantryAngles',oldRightLeafPoss',apertureInfoNew.beam(iNew).gantryAngle))';
                
                apertureInfoNew.beam(iNew).shape{phaseOld}(1).leftLeafPos_I = (interp1(oldGantryAngles',oldLeftLeafPoss',apertureInfoNew.propVMAT.beam(iNew).doseAngleBorders(1)))';
                apertureInfoNew.beam(iNew).shape{phaseOld}(1).rightLeafPos_I = (interp1(oldGantryAngles',oldRightLeafPoss',apertureInfoNew.propVMAT.beam(iNew).doseAngleBorders(1)))';
                
                apertureInfoNew.beam(iNew).shape{phaseOld}(1).leftLeafPos_F = (interp1(oldGantryAngles',oldLeftLeafPoss',apertureInfoNew.propVMAT.beam(iNew).doseAngleBorders(2)))';
                apertureInfoNew.beam(iNew).shape{phaseOld}(1).rightLeafPos_F = (interp1(oldGantryAngles',oldRightLeafPoss',apertureInfoNew.propVMAT.beam(iNew).doseAngleBorders(2)))';
            else
                apertureInfoNew.beam(iNew).shape{phaseOld}(1).leftLeafPos = apertureInfoOld.beam(iOld).shape{phaseOld}(1).leftLeafPos;
                apertureInfoNew.beam(iNew).shape{phaseOld}(1).rightLeafPos = apertureInfoOld.beam(iOld).shape{phaseOld}(1).rightLeafPos;
            end
            
            %all beams are now "optimized" to prevent their weights from being
            %overwritten
            %optAngleBorders becomes doseAngleBorders
            apertureInfoNew.beam(iNew).numOfShapes = 1;
            apertureInfoNew.propVMAT.beam(iNew).DAOBeam = true;
            %apertureInfoNew.beam(j).doseAngleOpt = stf(j).doseAngleOpt;
            apertureInfoNew.propVMAT.beam(iNew).DAOAngleBorders = stf(iNew).propVMAT.doseAngleBorders;
            apertureInfoNew.propVMAT.beam(iNew).DAOAngleBorderCentreDiff = stf(iNew).propVMAT.doseAngleBorderCentreDiff;
            apertureInfoNew.propVMAT.beam(iNew).DAOAngleBordersDiff = stf(iNew).propVMAT.doseAngleBordersDiff;
            apertureInfoNew.propVMAT.beam(iNew).timeFacCurr = apertureInfoNew.propVMAT.beam(iNew).doseAngleBordersDiff./apertureInfoNew.propVMAT.beam(iNew).DAOAngleBordersDiff; % = 1
            %apertureInfoNew.beam(j).timeFacPrev = stf(j).timeFacPrev;
            %apertureInfoNew.beam(j).timeFacNext = stf(j).timeFacNext;
            %apertureInfoNew.beam(j).IandFTimeInd = stf(j).IandFTimeInd;
            %apertureInfoNew.beam(j).IandFFac = stf(j).IandFFac;
            %apertureInfoNew.propVMAT.beam(j).timeFac = stf(j).propVMAT.timeFac;
            
            %{
        if isfield(recalc,'dijNew') && ~recalc.dijNew
            oldGantryAngles = [apertureInfoOld.beam.gantryAngle];
            diff = abs(apertureInfoNew.beam(j).gantryAngle-[apertureInfoOld.beam.gantryAngle]);
            minDiffInd = diff == min(diff);
            if isempty(stf(j).copyInd)
                recalc.pln.gantryAngles(j) = oldGantryAngles(minDiffInd);
            elseif stf(j).copyInd == 1
                recalc.pln.gantryAngles(j) = min(oldGantryAngles(minDiffInd));
            elseif stf(j).copyInd == 2
                recalc.pln.gantryAngles(j) = max(oldGantryAngles(minDiffInd));
            end
        end
            %}
            
            %{
        if apertureInfoOld.propVMAT.beam(i).FMOBeam
            apertureInfoNew.propVMAT.beam(j).FMOBeam = true;
            apertureInfoNew.propVMAT.beam(j).FMOAngleBorders = stf(j).propVMAT.FMOAngleBorders;
            apertureInfoNew.propVMAT.beam(j).FMOAngleBorderCentreDiff = stf(j).propVMAT.FMOAngleBorderCentreDiff;
            apertureInfoNew.propVMAT.beam(j).FMOAngleBordersDiff = stf(j).propVMAT.FMOAngleBordersDiff;
        else
            apertureInfoNew.propVMAT.beam(j).FMOBeam = false;
        end
            %}
            
            apertureInfoNew.apertureVector(shapeInd) = apertureInfoNew.beam(iNew).shape{phaseOld}(1).weight;
            shapeInd = shapeInd+1;
        end
        
    end
end

% do 4D stuff
if isfield(recalc,'numPhases')
    
    % update number of phases
    apertureInfoNew.numPhases = recalc.numPhases;
    
    % let all phases of the new plan library have the same set of apertures as
    % phase 1 in the old library
    for i = 1:numel(apertureInfoNew.beam)
        for phase = 1:apertureInfoNew.numPhases
            apertureInfoNew.beam(i).shape{phase} = apertureInfoNew.beam(i).shape{1};
        end
    end
    apertureInfoNew.motionModel.probPhase = recalc.pln.propOpt.prop4D.motionModel.probPhase;
end

apertureInfoNew.totalNumOfShapes        = sum([apertureInfoNew.beam.numOfShapes]);
apertureInfoNew.totalNumOfLeafPairs     = sum([apertureInfoNew.beam.numOfShapes]*[apertureInfoNew.beam.numOfActiveLeafPairs]');
apertureInfoNew.doseTotalNumOfLeafPairs = sum([apertureInfoNew.beam(:).numOfActiveLeafPairs]);
apertureInfoNew.totalNumOfOptBixels     = apertureInfoNew.totalNumOfBixels;

%recalc apertureVector
[apertureInfoNew.apertureVector, apertureInfoNew.mappingMx, apertureInfoNew.limMx] = matRad_daoApertureInfo2Vec(apertureInfoNew);

% recalc bixel weights
apertureInfoNew.propVMAT.continuousAperture = recalc.continuousAperture;
apertureInfoNew =  matRad_daoVec2ApertureInfo_VMATrecalcDynamic(apertureInfoNew,apertureInfoNew.apertureInfo.apertureVector);


recalc.apertureInfo = apertureInfoNew;
recalc.stf = stf;