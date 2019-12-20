function [apertureInfoVec, mappingMx, limMx] = matRad_daoApertureInfo2Vec(apertureInfo)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to generate a vector respresentation of the aperture
% weights and shapes and (optional) some meta information needed during
% optimization
%
% call
%   [apertureInfoVec, mappingMx, limMx] = matRad_daoApertureInfo2Vec(apertureInfo)
%
% input
%   apertureInfo:    aperture weight and shape info struct
%
% output
%   apertureInfoVec: vector respresentation of the apertue weights and shapes
%   mappingMx:       mapping of vector components to beams, shapes and leaves
%   limMx:           bounds on vector components, i.e., minimum and maximum
%                    aperture weights (0/inf) and leav positions (custom)
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

% function to create a single vector for the direct aperature optimization
% first: aperature weights (all phases)
% second: left leaf positions (all phases)
% third: right leaf positions (all phases)
% fourth (VMAT only): times between successive DAO gantry angles (all phases)

% initializing variables

vecLength = (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases; % each phase gets a weight and leaf pair

if apertureInfo.runVMAT && ~apertureInfo.propVMAT.fixedGantrySpeed
    vecLength = vecLength+apertureInfo.totalNumOfShapes; %Extra set of (apertureInfo.totalNumOfShapes) number of elements, allowing arc sector times to be optimized
end

apertureInfoVec = NaN * ones(vecLength,1);

offset = 0;

%% 1. aperture weights
for phase = 1:apertureInfo.numPhases
    for i = 1:size(apertureInfo.beam,2)
        for j = 1:apertureInfo.beam(i).numOfShapes
            
            %In VMAT, this weight is ~ "spread" over unoptimized beams (really, we assume constant dose rate DAO arc sector)
            apertureInfoVec(offset+j) = apertureInfo.beam(i).shape{phase}(j).jacobiScale*apertureInfo.beam(i).shape{phase}(j).weight;
            
        end
        offset = offset + apertureInfo.beam(i).numOfShapes;
    end
end

% 2. left and right leaf positions
%% fill the vector for all shapes of all beams
for phase = 1:apertureInfo.numPhases
    for i = 1:size(apertureInfo.beam,2)
        for j = 1:apertureInfo.beam(i).numOfShapes
            
            if ~apertureInfo.propVMAT.continuousAperture
                
                apertureInfoVec(offset+[1:apertureInfo.beam(i).numOfActiveLeafPairs]) = apertureInfo.beam(i).shape{phase}(j).leftLeafPos;
                apertureInfoVec(offset+[1:apertureInfo.beam(i).numOfActiveLeafPairs]+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases) = apertureInfo.beam(i).shape{phase}(j).rightLeafPos;
                
                offset = offset + apertureInfo.beam(i).numOfActiveLeafPairs;
            else
                % the leaf positions in the vector are defined at the
                % borders of the DAO arc, whereas they are defined at the
                % borders of the fluence arc in the apertureInfo struct 
                % (leftLeafPos_I, _F, etc.)
                % use fracFromLastDAOI_leafI etc. to convert between the
                % two
                % NOTE: here, since we are strictly using DAO beams,
                % lastDAOIndex = nextDAOIndex
                % so we can safely mix up the two
                
                % leafPos_I = fracFromLastDAOI_leafI.*leafPos_vecI+fracFromLastDAOF_leafI.*leafPos_vecF
                % leafPos_F = fracFromNextDAOI_leafF.*leafPos_vecI+fracFromNextDAOF_leafF.*leafPos_vecF
                
                fracToDAOI_leafI = apertureInfo.propVMAT.beam(i).fracFromNextDAOF_leafF ...
                    ./(apertureInfo.propVMAT.beam(i).fracFromLastDAOI_leafI.*apertureInfo.propVMAT.beam(i).fracFromNextDAOF_leafF-apertureInfo.propVMAT.beam(i).fracFromLastDAOF_leafI.*apertureInfo.propVMAT.beam(i).fracFromNextDAOI_leafF);
                fracToDAOI_leafF = -apertureInfo.propVMAT.beam(i).fracFromLastDAOF_leafI ...
                    ./(apertureInfo.propVMAT.beam(i).fracFromLastDAOI_leafI.*apertureInfo.propVMAT.beam(i).fracFromNextDAOF_leafF-apertureInfo.propVMAT.beam(i).fracFromLastDAOF_leafI.*apertureInfo.propVMAT.beam(i).fracFromNextDAOI_leafF);
                fracToDAOF_leafI = -apertureInfo.propVMAT.beam(i).fracFromNextDAOI_leafF ...
                    ./(apertureInfo.propVMAT.beam(i).fracFromLastDAOI_leafI.*apertureInfo.propVMAT.beam(i).fracFromNextDAOF_leafF-apertureInfo.propVMAT.beam(i).fracFromLastDAOF_leafI.*apertureInfo.propVMAT.beam(i).fracFromNextDAOI_leafF);
                fracToDAOF_leafF = apertureInfo.propVMAT.beam(i).fracFromLastDAOI_leafI ...
                    ./(apertureInfo.propVMAT.beam(i).fracFromLastDAOI_leafI.*apertureInfo.propVMAT.beam(i).fracFromNextDAOF_leafF-apertureInfo.propVMAT.beam(i).fracFromLastDAOF_leafI.*apertureInfo.propVMAT.beam(i).fracFromNextDAOI_leafF);
                
                if apertureInfo.propVMAT.beam(i).doseAngleDAO(1)
                    
                    % right then left
                    % ensure that the leaf positions are within bounds
                    leftLeafPos_vecI    = fracToDAOI_leafI.*apertureInfo.beam(i).shape{phase}(j).leftLeafPos_I+fracToDAOI_leafF.*apertureInfo.beam(i).shape{phase}(j).leftLeafPos_F;
                    rightLeafPos_vecI   = fracToDAOI_leafI.*apertureInfo.beam(i).shape{phase}(j).rightLeafPos_I+fracToDAOI_leafF.*apertureInfo.beam(i).shape{phase}(j).rightLeafPos_F;
                    
                    apertureInfoVec(offset+[1:apertureInfo.beam(i).numOfActiveLeafPairs]) = min(max(leftLeafPos_vecI,apertureInfo.beam(i).lim_l),apertureInfo.beam(i).lim_r);
                    apertureInfoVec(offset+[1:apertureInfo.beam(i).numOfActiveLeafPairs]+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases) = min(max(rightLeafPos_vecI,apertureInfo.beam(i).lim_l),apertureInfo.beam(i).lim_r);
                    
                    offset = offset + apertureInfo.beam(i).numOfActiveLeafPairs;
                end
                
                if apertureInfo.propVMAT.beam(i).doseAngleDAO(2)
                    
                    % right then left
                    % ensure that the leaf positions are within bounds
                    leftLeafPos_vecF    = fracToDAOF_leafI.*apertureInfo.beam(i).shape{phase}(j).leftLeafPos_I+fracToDAOF_leafF.*apertureInfo.beam(i).shape{phase}(j).leftLeafPos_F;
                    rightLeafPos_vecF   = fracToDAOF_leafI.*apertureInfo.beam(i).shape{phase}(j).rightLeafPos_I+fracToDAOF_leafF.*apertureInfo.beam(i).shape{phase}(j).rightLeafPos_F;
                    
                    apertureInfoVec(offset+[1:apertureInfo.beam(i).numOfActiveLeafPairs]) = min(max(leftLeafPos_vecF,apertureInfo.beam(i).lim_l),apertureInfo.beam(i).lim_r);
                    apertureInfoVec(offset+[1:apertureInfo.beam(i).numOfActiveLeafPairs]+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases) = min(max(rightLeafPos_vecF,apertureInfo.beam(i).lim_l),apertureInfo.beam(i).lim_r);
                    
                    offset = offset + apertureInfo.beam(i).numOfActiveLeafPairs;
                end
            end
            
        end
    end
end
%% 3. time of arc sector/beam
if apertureInfo.runVMAT && ~apertureInfo.propVMAT.fixedGantrySpeed
    offset = offset + apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
    
    %this gives a vector of the arc lengths belonging to each optimized CP
    %unique gets rid of double-counted angles (which is every interior
    %angle)
    
    optInd = [apertureInfo.propVMAT.beam.DAOBeam];
    optAngleLengths = [apertureInfo.propVMAT.beam(optInd).DAOAngleBordersDiff];
    optGantryRot = [apertureInfo.beam(optInd).gantryRot];
    apertureInfoVec((offset+1):end) = optAngleLengths./optGantryRot; %entries are the times until the next opt gantry angle is reached
    
end

%% 4. create additional information for later use
if nargout > 1
    
    mappingMx = NaN * ones(vecLength,5);
    limMx = NaN * ones(vecLength,2);
    
    limMx(1:(apertureInfo.totalNumOfShapes*apertureInfo.numPhases),:) = ones((apertureInfo.totalNumOfShapes*apertureInfo.numPhases),1)*[0 inf];
    
    counter = 1;
    for phase = 1:apertureInfo.numPhases
        for i = 1:numel(apertureInfo.beam)
            for j = 1:apertureInfo.beam(i).numOfShapes
                mappingMx(counter,1) = i;
                if apertureInfo.runVMAT && phase == 1 && ~apertureInfo.propVMAT.fixedGantrySpeed
                    fileName = apertureInfo.propVMAT.machineConstraintFile;
                    try
                        load(fileName,'machine');
                    catch
                        error(['Could not find the following machine file: ' fileName ]);
                    end
                    
                    timeLimL = diff(apertureInfo.propVMAT.beam(i).DAOAngleBorders)/machine.constraints.gantryRotationSpeed(2); %Minimum time interval between two optimized beams/gantry angles
                    timeLimU = diff(apertureInfo.propVMAT.beam(i).DAOAngleBorders)/machine.constraints.gantryRotationSpeed(1); %Maximum time interval between two optimized beams/gantry angles
                    
                    mappingMx(counter+(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases,1) = i;
                    limMx(counter+(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases,:) = [timeLimL timeLimU];
                end
                counter = counter + 1;
            end
        end
    end
    
    shapeOffset = 0;
    for phase = 1:apertureInfo.numPhases
        for i = 1:numel(apertureInfo.beam)
            for j = 1:apertureInfo.beam(i).numOfShapes
                for k = 1:apertureInfo.beam(i).numOfActiveLeafPairs
                    mappingMx(counter,1) = i;
                    mappingMx(counter,2) = j + shapeOffset; % store global shape number for grad calc
                    mappingMx(counter,3) = j; % store local shape number
                    mappingMx(counter,4) = k; % store local leaf number
                    mappingMx(counter,5) = phase; % store phase number
                    
                    limMx(counter,1)     = apertureInfo.beam(i).lim_l(k);
                    limMx(counter,2)     = apertureInfo.beam(i).lim_r(k);
                    counter = counter + 1;
                    
                    if apertureInfo.runVMAT && apertureInfo.propVMAT.continuousAperture && nnz(apertureInfo.propVMAT.beam(i).doseAngleDAO) == 2
                        %redo for initial and final leaf positions
                        %might have to revisit this after looking at gradient,
                        %esp. mappingMx(counter,2)
                        %only an issue for non-interpolated deliveries
                        mappingMx(counter,1) = i;
                        mappingMx(counter,2) = j + shapeOffset; % store global shape number for grad calc
                        mappingMx(counter,3) = j; % store local shape number
                        mappingMx(counter,4) = k; % store local leaf number
                        mappingMx(counter,5) = phase; % store phase number
                        
                        limMx(counter,1)     = apertureInfo.beam(i).lim_l(k);
                        limMx(counter,2)     = apertureInfo.beam(i).lim_r(k);
                        counter = counter + 1;
                    end
                end
            end
            shapeOffset = shapeOffset + apertureInfo.beam(i).numOfShapes;
        end
    end
    
end

mappingMx(counter:(apertureInfo.numPhases*(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)),:) = mappingMx(apertureInfo.numPhases*apertureInfo.totalNumOfShapes+1:counter-1,:);
limMx(counter:(apertureInfo.numPhases*(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)),:)     = limMx(apertureInfo.numPhases*apertureInfo.totalNumOfShapes+1:counter-1,:);

end


