function apertureInfo = matRad_sequencing2ApertureInfo(sequencing,stf,pln)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to generate a shape info struct based on the result of
% multileaf collimator sequencing
%
% call
%   apertureInfo = matRad_sequencing2ApertureInfo(Sequencing,stf)
%
% input
%   sequencing: matRad sequencing result struct
%   stf:        matRad steering information struct
%   pln:        matRad pln struct
%
% output
%   apertureInfo: matRad aperture weight and shape info struct
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


% MLC parameters:
bixelWidth = pln.propStf.bixelWidth; % [mm]
numOfMLCLeafPairs = 80;
%     define central leaf pair (here we want the 0mm position to be in the
%     center of a leaf pair (e.g. leaf 41 stretches from -2.5mm to 2.5mm
%     for a bixel/leafWidth of 5mm and 81 leaf pairs)
centralLeafPair = ceil(numOfMLCLeafPairs/2);

% initializing variables
bixelIndOffset = 0; % used for creation of bixel index maps
totalNumOfBixels = sum([stf(:).totalNumOfBixels]);
totalNumOfShapes = sum([sequencing.beam.numOfShapes]);
if pln.propOpt.VMAToptions.continuousAperture
    totalNumOfShapes = totalNumOfShapes-1;
end
weightOffset = 1;
vectorOffset = totalNumOfShapes + 1; % used for bookkeeping in the vector for optimization
bixOffset = 1; %used for gradient calculations

if pln.propOpt.runVMAT
    totalNumOfOptBixels = 0;
    totalNumOfLeafPairs = 0;
    
    apertureInfo.propVMAT.jacobT = zeros(totalNumOfShapes,numel(sequencing.beam));
    
    % preallocate propVMAT.beam (metadata)
    stf(stf(1).propVMAT.beamParentIndex).propVMAT.leafDir = 1;
    apertureInfo.propVMAT.beam(size(stf,2)) = stf(stf(1).propVMAT.beamParentIndex).propVMAT;
end

apertureInfo.jacobiScale = ones(totalNumOfShapes.*sequencing.numPhases,1);

% loop over all beams
for i = 1:size(stf,2)
    
    %% 1. read stf and create maps (Ray & Bixelindex)
    
    % get x- and z-coordinates of bixels
    rayPos_bev = reshape([stf(i).ray.rayPos_bev],3,[]);
    X = rayPos_bev(1,:)';
    Z = rayPos_bev(3,:)';
    
    % create ray-map
    maxX = max(X); minX = min(X);
    maxZ = max(Z); minZ = min(Z);
    
    dimX = (maxX-minX)/bixelWidth + 1;
    dimZ = (maxZ-minZ)/bixelWidth + 1;
    
    rayMap = zeros(dimZ,dimX);
    
    % get indices for x and z positions
    xPos = (X-minX)/bixelWidth + 1;
    zPos = (Z-minZ)/bixelWidth + 1;
    
    % get indices in the ray-map
    indInRay = zPos + (xPos-1)*dimZ;
    
    % fill ray-map
    rayMap(indInRay) = 1;
    
    if stf(i).propVMAT.doseBeam
        % for fluence-only beams, we will determine bixel index maps later
        
        % create map of bixel indices
        bixelIndMap = NaN * ones(dimZ,dimX);
        bixelIndMap(indInRay) = (1:stf(i).numOfRays) + bixelIndOffset;
        bixelIndOffset = bixelIndOffset + stf(i).numOfRays;
    end
    
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
    
    % get leaf positions for all shapes
    % leaf positions can be extracted from the shapes created in Sequencing
    for m = 1:sequencing.beam(i).numOfShapes
        
        % loading shape from Sequencing result
        shapeMap = sequencing.beam(i).shapes(:,:,m);
        % get left and right leaf indices from shapemap
        % initializing limits
        leftLeafPos = NaN * ones(dimZ,1);
        rightLeafPos = NaN * ones(dimZ,1);
        % looping over leaf pairs
        for l = 1:dimZ
            leftLeafPosInd  = find(shapeMap(l,:),1,'first');
            rightLeafPosInd = find(shapeMap(l,:),1,'last');
            
            if isempty(leftLeafPosInd) && isempty(rightLeafPosInd) % if no bixel is open, use limits from Ray positions
                leftLeafPos(l) = (lim_l(l)+lim_r(l))/2;
                rightLeafPos(l) = leftLeafPos(l);
            else
                % the physical position [mm] can be calculated from the indices
                leftLeafPos(l) = (leftLeafPosInd-1)*bixelWidth...
                    + minX - 1/2*bixelWidth;
                rightLeafPos(l) = (rightLeafPosInd-1)*bixelWidth...
                    + minX + 1/2*bixelWidth;
                
                %Can happen in some cases in SW trajectory sampling
                if leftLeafPos(l) < lim_l(l)
                    leftLeafPos(l) = lim_l(l);
                end
                if rightLeafPos(l) > lim_r(l)
                    rightLeafPos(l) = lim_r(l);
                end
                
            end
        end
        
        if pln.propOpt.runVMAT && pln.propOpt.VMAToptions.continuousAperture
            
            if sequencing.beam(i).numOfShapes == 1
                % this is not the first DAO beam
                
                % this data is saved to _F_DAO (to be eliminated in
                % daoApertureInfo2Vec)
                apertureInfo.beam(i).shape{1}.leftLeafPos_F_DAO = leftLeafPos;
                apertureInfo.beam(i).shape{1}.rightLeafPos_F_DAO = rightLeafPos;
                apertureInfo.beam(i).shape{1}.weight = sequencing.beam(i).shapesWeight;
                apertureInfo.beam(i).shape{1}.shapeMap = shapeMap;
            elseif sequencing.beam(i).numOfShapes == 2
                % this is the first DAO beam
                
                % change numOfShapes to 1
                sequencing.beam(i).numOfShapes = 1;
                
                % m = 1 is saved to _I_DAO, m = 2 to _F_DAO (both are to be
                % eliminated in daoApertureInfo2Vec)
                if m == 1
                    apertureInfo.beam(i).shape{1}.leftLeafPos_I_DAO = leftLeafPos;
                    apertureInfo.beam(i).shape{1}.rightLeafPos_I_DAO = rightLeafPos;
                    apertureInfo.beam(i).shape{1}.weight = sequencing.beam(i).shapesWeight;
                    apertureInfo.beam(i).shape{1}.shapeMap = shapeMap;
                elseif m == 2
                    apertureInfo.beam(i).shape{1}.leftLeafPos_F_DAO = leftLeafPos;
                    apertureInfo.beam(i).shape{1}.rightLeafPos_F_DAO = rightLeafPos;
                    apertureInfo.beam(i).shape{1}.weight = sequencing.beam(i).shapesWeight;
                    apertureInfo.beam(i).shape{1}.shapeMap = shapeMap;
                end
            end
        else
            % save data for each shape of this beam
            apertureInfo.beam(i).shape{1}(m).leftLeafPos = leftLeafPos;
            apertureInfo.beam(i).shape{1}(m).rightLeafPos = rightLeafPos;
            apertureInfo.beam(i).shape{1}(m).weight = sequencing.beam(i).shapesWeight(m);
            apertureInfo.beam(i).shape{1}(m).shapeMap = shapeMap;
        end
        
        if pln.propOpt.runVMAT
            apertureInfo.beam(i).shape{1}.MURate = sequencing.beam(i).MURate;
            if m == 2
                continue
            end
        end
        
        apertureInfo.beam(i).shape{1}(m).jacobiScale = 1;
        
        if pln.propOpt.VMAToptions.continuousAperture
            apertureInfo.beam(i).shape{1}(m).vectorOffset = [vectorOffset vectorOffset+dimZ];
            
            % update index for bookkeeping
            vectorOffset = vectorOffset + dimZ*nnz(stf(i).propVMAT.doseAngleDAO);
        else
            apertureInfo.beam(i).shape{1}(m).vectorOffset = [vectorOffset vectorOffset];
            
            % update index for bookkeeping
            vectorOffset = vectorOffset + dimZ;
        end
        apertureInfo.beam(i).shape{1}(m).weightOffset = weightOffset;
        weightOffset = weightOffset+1;
    end
    
    % z-coordinates of active leaf pairs
    % get z-coordinates from bixel positions
    leafPairPos = unique(Z);
    
    % find upmost and downmost leaf pair
    topLeafPairPos = maxZ;
    bottomLeafPairPos = minZ;
    
    topLeafPair = centralLeafPair - topLeafPairPos/bixelWidth;
    bottomLeafPair = centralLeafPair - bottomLeafPairPos/bixelWidth;
    
    % create bool map of active leaf pairs
    isActiveLeafPair = zeros(numOfMLCLeafPairs,1);
    isActiveLeafPair(topLeafPair:bottomLeafPair) = 1;
    
    % create MLC window
    % getting the dimensions of the MLC in order to be able to plot the
    % shapes using physical coordinates
    MLCWindow = [minX-bixelWidth/2 maxX+bixelWidth/2 ...
        minZ-bixelWidth/2 maxZ+bixelWidth/2];
    
    % save data for each beam
    apertureInfo.beam(i).numOfShapes = sequencing.beam(i).numOfShapes;
    apertureInfo.beam(i).numOfActiveLeafPairs = dimZ;
    apertureInfo.beam(i).leafPairPos = leafPairPos;
    apertureInfo.beam(i).isActiveLeafPair = isActiveLeafPair;
    apertureInfo.beam(i).centralLeafPair = centralLeafPair;
    apertureInfo.beam(i).lim_l = lim_l;
    apertureInfo.beam(i).lim_r = lim_r;
    apertureInfo.beam(i).posOfCornerBixel = posOfCornerBixel;
    apertureInfo.beam(i).MLCWindow = MLCWindow;
    apertureInfo.beam(i).gantryAngle = stf(i).gantryAngle;
    if stf(i).propVMAT.doseBeam
        % for fluence-only beams, we will determine bixel index maps later
        apertureInfo.beam(i).bixelIndMap = bixelIndMap;
        apertureInfo.beam(i).numBixels = nnz(~isnan(bixelIndMap));
    end
    
    if pln.propOpt.runVMAT
        
        apertureInfo.beam(i).bixOffset = bixOffset;
        bixOffset = bixOffset+apertureInfo.beam(i).numOfActiveLeafPairs;
        
        stf(i).propVMAT.leafDir = sequencing.beam(i).leafDir;
        
        % put all propVMAT stuff from stf into apertureInfo
        apertureInfo.propVMAT.beam(i) = stf(i).propVMAT;
        
        apertureInfo.beam(i).gantryRot  = sequencing.beam(i).gantryRot;
        apertureInfo.beam(i).time       = apertureInfo.propVMAT.beam(i).fluAngleBordersDiff./apertureInfo.beam(i).gantryRot;
        
        if apertureInfo.propVMAT.beam(i).DAOBeam
            
            totalNumOfOptBixels = totalNumOfOptBixels+stf(i).totalNumOfBixels;
            totalNumOfLeafPairs = totalNumOfLeafPairs+apertureInfo.beam(i).numOfShapes*apertureInfo.beam(i).numOfActiveLeafPairs;
            
            apertureInfo.propVMAT.jacobT(stf(i).propVMAT.DAOIndex,i) = stf(i).propVMAT.timeFacCurr;
            
        else
            
            % gradient of the time in the current fluence calculation arc
            % wrt optimized DAO arc times
            % it is equal to the fractional portion of the fluence arc
            % spent in the last/next DAO arcs
            apertureInfo.propVMAT.jacobT(stf(stf(i).propVMAT.lastDAOIndex).propVMAT.DAOIndex,i) = apertureInfo.propVMAT.jacobT(stf(stf(i).propVMAT.lastDAOIndex).propVMAT.DAOIndex,i)+stf(stf(i).propVMAT.lastDAOIndex).propVMAT.timeFacCurr.*stf(i).propVMAT.fracFromLastDAO_gantryRot.*stf(i).propVMAT.fluAngleBordersDiff./stf(stf(i).propVMAT.lastDAOIndex).propVMAT.fluAngleBordersDiff;
            apertureInfo.propVMAT.jacobT(stf(stf(i).propVMAT.nextDAOIndex).propVMAT.DAOIndex,i) = apertureInfo.propVMAT.jacobT(stf(stf(i).propVMAT.nextDAOIndex).propVMAT.DAOIndex,i)+stf(stf(i).propVMAT.nextDAOIndex).propVMAT.timeFacCurr.*stf(i).propVMAT.fracFromNextDAO_gantryRot.*stf(i).propVMAT.fluAngleBordersDiff./stf(stf(i).propVMAT.nextDAOIndex).propVMAT.fluAngleBordersDiff;
        end
    end
end

% save global data
apertureInfo.runVMAT                = pln.propOpt.runVMAT;
apertureInfo.preconditioner         = pln.propOpt.preconditioner;
apertureInfo.run4D                  = pln.propOpt.run4D;
apertureInfo.varOpt                 = pln.propOpt.varOpt;
apertureInfo.numPhases              = sequencing.numPhases;
apertureInfo.bixelWidth             = bixelWidth;
apertureInfo.numOfMLCLeafPairs      = numOfMLCLeafPairs;
apertureInfo.totalNumOfBixels       = totalNumOfBixels;
apertureInfo.totalNumOfOptBixels    = totalNumOfOptBixels;
apertureInfo.totalNumOfShapes       = sum([apertureInfo.beam.numOfShapes]);

if isfield(sequencing,'weightToMU')
    apertureInfo.weightToMU = sequencing.weightToMU;
end

% save more metadata
apertureInfo = matRad_apertureInfoMeta(apertureInfo,pln,stf);

% create vectors for optimization
[apertureInfo.apertureVector, apertureInfo.mappingMx, apertureInfo.limMx] = matRad_daoApertureInfo2Vec(apertureInfo);

end

