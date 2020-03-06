function stf = matRad_StfVMATPost(stf,pln,masterRayPosBEV,masterTargetPointBEV,SAD,machine)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad steering information post-processing for VMAT
%
% call
%   stf = matRad_StfVMATPost(stf,pln)
%
% input
%   stf:        matRad steering information struct
%   pln:        matRad plan meta information struct
%
% output
%   stf:        matRad steering information struct
%
% References
%   -
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

%After all steering file information is completed, loop over
%initialized gantry angles.  All children and subchildren of these angles should
%have ray positions given by the union of their own ray positions and
%the ray positions of the parent transformed to their gantry angle.
%This is so that: (1) the target is still totally in the FOV of each
%angle; and (2) the parent can give segments to the children during
%initial segmentation and DAO.

fprintf('matRad: VMAT post-processing (1/2)... ');

numDAO = 1;
DAOAngleBorders = zeros(2*numel(pln.propStf.DAOGantryAngles),1);
offset = 1;
timeFacIndOffset = 1;
firstDAO = true;
firstFMO = true;

% initialize values for last/next DAO/dose indices
lastDAOIndex    = find(abs(pln.propStf.fluGantryAngles - pln.propStf.DAOGantryAngles(1)) < 1e-8);
nextDAOIndex    = find(abs(pln.propStf.fluGantryAngles - pln.propStf.DAOGantryAngles(2)) < 1e-8);
lastDoseIndex   = find(abs(pln.propStf.fluGantryAngles - pln.propStf.gantryAngles(1)) < 1e-8);
nextDoseIndex   = find(abs(pln.propStf.fluGantryAngles - pln.propStf.gantryAngles(2)) < 1e-8);

for i = 1:length(stf)
    
    % Determine which FMO beam the current beam belongs to
    [~,stf(i).propVMAT.beamParentFMOIndex] = min(abs(pln.propStf.FMOGantryAngles-pln.propStf.fluGantryAngles(i)));
    stf(i).propVMAT.beamParentGantryAngle = pln.propStf.FMOGantryAngles(stf(i).propVMAT.beamParentFMOIndex);
    stf(i).propVMAT.beamParentIndex = find(abs(pln.propStf.fluGantryAngles - stf(i).propVMAT.beamParentGantryAngle) < 1e-6);
    
    % Indicate if this beam is to be included in DOA/FMO/dosecalc or not.
    % NB: not all beams are considered in dose calc for objective function
    % in DAO: some are only used for fluence calculation
    stf(i).propVMAT.FMOBeam     = any(abs(pln.propStf.FMOGantryAngles - pln.propStf.fluGantryAngles(i)) < 1e-6);
    stf(i).propVMAT.DAOBeam     = any(abs(pln.propStf.DAOGantryAngles - pln.propStf.fluGantryAngles(i)) < 1e-6);
    stf(i).propVMAT.doseBeam    = any(abs(pln.propStf.gantryAngles - pln.propStf.fluGantryAngles(i)) < 1e-6);
    
    if stf(i).propVMAT.DAOBeam
        if firstDAO
            stf(i).propVMAT.firstDAO = true;
            firstDAO = false;
        else
            stf(i).propVMAT.firstDAO = false;
        end
    else
        % dummy variables for DAO stuff
        stf(i).propVMAT.firstDAO = [];
    end
    
    if stf(i).propVMAT.FMOBeam
        if firstFMO
            stf(i).propVMAT.firstFMO = true;
            firstFMO = false;
        else
            stf(i).propVMAT.firstFMO = false;
        end
    else
        % dummy variables for FMO stuff
        stf(i).propVMAT.numOfBeamSubChildren        = [];
        stf(i).propVMAT.beamSubChildrenGantryAngles = [];
        stf(i).propVMAT.beamSubChildrenIndex        = [];
        stf(i).propVMAT.numOfBeamChildren           = [];
        stf(i).propVMAT.beamChildrenGantryAngles    = [];
        stf(i).propVMAT.beamChildrenIndex           = [];
        stf(i).propVMAT.firstFMO                    = [];
    end
    
    %% Determine different angle borders
    
    % fluAngleBorders are the angular borders over which fluence is deposited
    if i == 1
        
        stf(i).propVMAT.fluAngleBorders = [pln.propOpt.VMAToptions.startingAngle (pln.propStf.fluGantryAngles(i+1)+pln.propStf.fluGantryAngles(i))/2];
    elseif i == length(pln.propStf.fluGantryAngles)
        
        stf(i).propVMAT.fluAngleBorders = [(pln.propStf.fluGantryAngles(i-1)+pln.propStf.fluGantryAngles(i))/2 pln.propOpt.VMAToptions.finishingAngle];
    else
        
        stf(i).propVMAT.fluAngleBorders = ([pln.propStf.fluGantryAngles(i-1) pln.propStf.fluGantryAngles(i+1)]+pln.propStf.fluGantryAngles(i))/2;
    end
    
    stf(i).propVMAT.fluAngleBorderCentreDiff = [stf(i).gantryAngle-stf(i).propVMAT.fluAngleBorders(1) stf(i).propVMAT.fluAngleBorders(2)-stf(i).gantryAngle];
    stf(i).propVMAT.fluAngleBordersDiff = sum(stf(i).propVMAT.fluAngleBorderCentreDiff);
    
    if stf(i).propVMAT.doseBeam
        % doseAngleBorders are the angular borders over which dose is deposited
        doseIndex = find(abs(pln.propStf.gantryAngles - pln.propStf.fluGantryAngles(i)) < 1e-8);
        
        if doseIndex == 1
            stf(i).propVMAT.doseAngleBorders = [pln.propOpt.VMAToptions.startingAngle (pln.propStf.gantryAngles(doseIndex+1)+pln.propStf.gantryAngles(doseIndex))/2];
            
            lastDoseIndex = i;
            nextDoseIndex = find(abs(pln.propStf.fluGantryAngles - pln.propStf.gantryAngles(doseIndex+1)) < 1e-8);
        elseif doseIndex == length(pln.propStf.gantryAngles)
            stf(i).propVMAT.doseAngleBorders = [(pln.propStf.gantryAngles(doseIndex-1)+pln.propStf.gantryAngles(doseIndex))/2 pln.propOpt.VMAToptions.finishingAngle];
            
            lastDoseIndex = find(abs(pln.propStf.fluGantryAngles - pln.propStf.gantryAngles(doseIndex-1)) < 1e-8);
            nextDoseIndex = i;
        else
            stf(i).propVMAT.doseAngleBorders = ([pln.propStf.gantryAngles(doseIndex-1) pln.propStf.gantryAngles(doseIndex+1)]+pln.propStf.gantryAngles(doseIndex))/2;
            
            lastDoseIndex = i;
            nextDoseIndex = find(abs(pln.propStf.fluGantryAngles - pln.propStf.gantryAngles(doseIndex+1)) < 1e-8);
        end
        
        stf(i).propVMAT.doseAngleBorderCentreDiff = [stf(i).gantryAngle-stf(i).propVMAT.doseAngleBorders(1) stf(i).propVMAT.doseAngleBorders(2)-stf(i).gantryAngle];
        stf(i).propVMAT.doseAngleBordersDiff = sum(stf(i).propVMAT.doseAngleBorderCentreDiff);
        
        
        if stf(i).propVMAT.FMOBeam
            % FMOAngleBorders are the angular borders over which an optimized
            % control point has influence
            FMOIndex = find(abs(pln.propStf.FMOGantryAngles - pln.propStf.fluGantryAngles(i)) < 1e-8);
            
            if FMOIndex == 1
                
                stf(i).propVMAT.FMOAngleBorders = [pln.propOpt.VMAToptions.startingAngle (pln.propStf.FMOGantryAngles(FMOIndex+1)+pln.propStf.FMOGantryAngles(FMOIndex))/2];
            elseif FMOIndex == length(pln.propStf.FMOGantryAngles)
                
                stf(i).propVMAT.FMOAngleBorders = [(pln.propStf.FMOGantryAngles(FMOIndex-1)+pln.propStf.FMOGantryAngles(FMOIndex))/2 pln.propOpt.VMAToptions.finishingAngle];
            else
                
                stf(i).propVMAT.FMOAngleBorders = ([pln.propStf.FMOGantryAngles(FMOIndex-1) pln.propStf.FMOGantryAngles(FMOIndex+1)]+pln.propStf.FMOGantryAngles(FMOIndex))/2;
            end
            stf(i).propVMAT.FMOAngleBorderCentreDiff = [stf(i).gantryAngle-stf(i).propVMAT.FMOAngleBorders(1) stf(i).propVMAT.FMOAngleBorders(2)-stf(i).gantryAngle];
            stf(i).propVMAT.FMOAngleBordersDiff = sum(stf(i).propVMAT.FMOAngleBorderCentreDiff);
        else
            % dummy variables for FMOAngleBorders
            stf(i).propVMAT.FMOAngleBorders             = [];
            stf(i).propVMAT.FMOAngleBorderCentreDiff    = [];
            stf(i).propVMAT.FMOAngleBordersDiff         = [];
        end
    else
        % dummy variables for doseAngleBorders
        stf(i).propVMAT.doseAngleBorders            = [];
        stf(i).propVMAT.doseAngleBorderCentreDiff   = [];
        stf(i).propVMAT.doseAngleBordersDiff        = [];
        
        % dummy variables for FMOAngleBorders
        stf(i).propVMAT.FMOAngleBorders             = [];
        stf(i).propVMAT.FMOAngleBorderCentreDiff    = [];
        stf(i).propVMAT.FMOAngleBordersDiff         = [];
    end
    
    %Assign beam to its Parent, either as child (optimized) or subchild
    %(interpolated)
    if stf(i).propVMAT.DAOBeam
        
        if ~isfield(stf(stf(i).propVMAT.beamParentIndex).propVMAT,'beamChildrenGantryAngles') || isempty(stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamChildrenGantryAngles)
            stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamChildren = 0;
            stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamChildrenGantryAngles = nan(1000,1);
            stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamChildrenIndex = nan(1000,1);
        end
        
        stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamChildren = stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamChildren+1;
        stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamChildrenGantryAngles(stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamChildren) = pln.propStf.fluGantryAngles(i);
        stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamChildrenIndex(stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamChildren) = i;
        
        %optAngleBorders are the angular borders over which an optimized control point
        %has influence
        DAOIndex = find(abs(pln.propStf.DAOGantryAngles - pln.propStf.fluGantryAngles(i)) < 1e-8);
        
        if DAOIndex == 1
            stf(i).propVMAT.DAOAngleBorders = [pln.propOpt.VMAToptions.startingAngle (pln.propStf.DAOGantryAngles(DAOIndex+1)+pln.propStf.DAOGantryAngles(DAOIndex))/2];
            
            lastDAOIndex = i;
            nextDAOIndex = find(abs(pln.propStf.fluGantryAngles - pln.propStf.DAOGantryAngles(DAOIndex+1)) < 1e-8);
        elseif DAOIndex == length(pln.propStf.DAOGantryAngles)
            stf(i).propVMAT.DAOAngleBorders = [(pln.propStf.DAOGantryAngles(DAOIndex-1)+pln.propStf.DAOGantryAngles(DAOIndex))/2 pln.propOpt.VMAToptions.finishingAngle];
            
            lastDAOIndex = find(abs(pln.propStf.fluGantryAngles - pln.propStf.DAOGantryAngles(DAOIndex-1)) < 1e-8);
            nextDAOIndex = i;
        else
            stf(i).propVMAT.DAOAngleBorders = ([pln.propStf.DAOGantryAngles(DAOIndex-1) pln.propStf.DAOGantryAngles(DAOIndex+1)]+pln.propStf.DAOGantryAngles(DAOIndex))/2;
            
            lastDAOIndex = i;
            nextDAOIndex = find(abs(pln.propStf.fluGantryAngles - pln.propStf.DAOGantryAngles(DAOIndex+1)) < 1e-8);
        end
        
        % collect all DAO angle borders
        DAOAngleBorders((offset):(offset+1)) = stf(i).propVMAT.DAOAngleBorders;
        offset= offset+2;
        
        stf(i).propVMAT.DAOIndex = numDAO;
        numDAO = numDAO+1;
        
        stf(i).propVMAT.DAOAngleBorderCentreDiff = [stf(i).gantryAngle-stf(i).propVMAT.DAOAngleBorders(1) stf(i).propVMAT.DAOAngleBorders(2)-stf(i).gantryAngle];
        stf(i).propVMAT.DAOAngleBordersDiff = sum(stf(i).propVMAT.DAOAngleBorderCentreDiff);
        
        %This is the factor that relates the total time in the
        %optimized arc sector to the total time in the current dose
        %sector
        stf(i).propVMAT.timeFacCurr =  stf(i).propVMAT.fluAngleBordersDiff./stf(i).propVMAT.DAOAngleBordersDiff;
        
        if pln.propOpt.VMAToptions.continuousAperture
            %These are the factors that relate the total time in the
            %optimized arc sector to the total time in the previous and
            %next dose sectors
            stf(i).propVMAT.timeFac1 = zeros(1,3);
            
            stf(i).propVMAT.timeFac1(1) = (stf(i).propVMAT.DAOAngleBorderCentreDiff(1)-stf(i).propVMAT.fluAngleBorderCentreDiff(1))/stf(i).propVMAT.DAOAngleBordersDiff;
            stf(i).propVMAT.timeFac1(2) = stf(i).propVMAT.timeFacCurr;
            stf(i).propVMAT.timeFac1(3) = (stf(i).propVMAT.DAOAngleBorderCentreDiff(2)-stf(i).propVMAT.fluAngleBorderCentreDiff(2))/stf(i).propVMAT.DAOAngleBordersDiff;
            
            % keep entries with a non-0 timeFac
            delInd     = stf(i).propVMAT.timeFac1 == 0;
            
            % write timeFacInd
            stf(i).propVMAT.timeFacInd1          = [timeFacIndOffset-1 timeFacIndOffset timeFacIndOffset+1];
            stf(i).propVMAT.timeFacInd1(delInd)  = 0;
            
            % update offset
            if delInd(3)
                timeFacIndOffset = timeFacIndOffset+1;
            else
                timeFacIndOffset = timeFacIndOffset+2;
            end
            
        else
            %These are the factors that relate the total time in the
            %optimized arc sector to the total time in the previous and
            %next dose sectors
            stf(i).propVMAT.timeFac = zeros(1,2);
            
            stf(i).propVMAT.timeFac(1) = stf(i).propVMAT.DAOAngleBorderCentreDiff(1)/stf(i).propVMAT.DAOAngleBordersDiff;
            stf(i).propVMAT.timeFac(2) = stf(i).propVMAT.DAOAngleBorderCentreDiff(2)/stf(i).propVMAT.DAOAngleBordersDiff;
        end
        
    else
        if ~isfield(stf(stf(i).propVMAT.beamParentIndex).propVMAT,'beamSubChildrenGantryAngles') || isempty(stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamSubChildrenGantryAngles)
            stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamSubChildren = 0;
            stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamSubChildrenGantryAngles = nan(1000,1);
            stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamSubChildrenIndex = nan(1000,1);
        end
        
        stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamSubChildren = stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamSubChildren+1;
        stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamSubChildrenGantryAngles(stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamSubChildren) = pln.propStf.fluGantryAngles(i);
        stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamSubChildrenIndex(stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamSubChildren) = i;
        
        % dummy variables for DAOAngleBorders
        stf(i).propVMAT.DAOAngleBorders             = [];
        stf(i).propVMAT.DAOAngleBorderCentreDiff    = [];
        stf(i).propVMAT.DAOAngleBordersDiff         = [];
        
        % dummy variables for DAOIndex
        stf(i).propVMAT.DAOIndex = [];
        
        % dummy variables for timeFac
        stf(i).propVMAT.timeFacCurr = [];
        stf(i).propVMAT.timeFac1    = [];
        stf(i).propVMAT.timeFacInd1 = [];
    end
    
    % store last/next DAO/dose indices
    stf(i).propVMAT.lastDAOIndex    = lastDAOIndex;
    stf(i).propVMAT.nextDAOIndex    = nextDAOIndex;
    stf(i).propVMAT.lastDoseIndex   = lastDoseIndex;
    stf(i).propVMAT.nextDoseIndex   = nextDoseIndex;
    
    %% transformation of union of rays
    
    % only necessary for dose beams
    if stf(i).propVMAT.doseBeam
        
        stf(i).numOfRays = size(masterRayPosBEV,1);
        stf(i).numOfBixelsPerRay = ones(1,stf(i).numOfRays);
        stf(i).totalNumOfBixels = sum(stf(i).numOfBixelsPerRay);
        
        
        % source position in bev
        stf(i).sourcePoint_bev = [0 -SAD 0];
        
        % get (active) rotation matrix
        % transpose matrix because we are working with row vectors
        rotMat_vectors_T = transpose(matRad_getRotationMatrix(stf(i).gantryAngle,stf(i).couchAngle));
        
        
        stf(i).sourcePoint = stf(i).sourcePoint_bev*rotMat_vectors_T;
        
        % Save ray and target position in lps system.
        for j = 1:stf(i).numOfRays
            stf(i).ray(j).rayPos_bev = masterRayPosBEV(j,:);
            stf(i).ray(j).targetPoint_bev = masterTargetPointBEV(j,:);
            
            stf(i).ray(j).rayPos      = stf(i).ray(j).rayPos_bev*rotMat_vectors_T;
            stf(i).ray(j).targetPoint = stf(i).ray(j).targetPoint_bev*rotMat_vectors_T;
            if strcmp(pln.radiationMode,'photons')
                stf(i).ray(j).rayCorners_SCD = (repmat([0, machine.meta.SCD - SAD, 0],4,1)+ (machine.meta.SCD/SAD) * ...
                    [masterRayPosBEV(j,:) + [+stf(i).bixelWidth/2,0,+stf(i).bixelWidth/2];...
                    masterRayPosBEV(j,:) + [-stf(i).bixelWidth/2,0,+stf(i).bixelWidth/2];...
                    masterRayPosBEV(j,:) + [-stf(i).bixelWidth/2,0,-stf(i).bixelWidth/2];...
                    masterRayPosBEV(j,:) + [+stf(i).bixelWidth/2,0,-stf(i).bixelWidth/2]])*rotMat_vectors_T;
            end
        end
        
        % loop over all rays to determine meta information for each ray
        stf(i).totalNumOfBixels = sum(stf(i).numOfBixelsPerRay);
        stf(i).numOfBixelsPerRay = ones(1,stf(i).numOfRays);
        
        for j = stf(i).numOfRays:-1:1
            
            % find appropriate energies for particles
            if strcmp(stf(i).radiationMode,'photons')
                
                % book keeping for photons
                stf(i).ray(j).energy = machine.data.energy;
            else
                error('Error generating stf struct: invalid radiation modality for VMAT.');
            end
        end
    end
    
    matRad_progress(i,length(stf));
end

%% final cleanup and calculation of factors we couldn't calc before
fprintf('matRad: VMAT post-processing (2/2)... ');

for i = 1:length(stf)
    if stf(i).propVMAT.FMOBeam
        %remove NaNs from beamChildren and beamSubChildren
        if isfield(stf(i).propVMAT,'beamChildrenGantryAngles')
            stf(i).propVMAT.beamChildrenGantryAngles(isnan(stf(i).propVMAT.beamChildrenGantryAngles)) = [];
            stf(i).propVMAT.beamChildrenIndex(isnan(stf(i).propVMAT.beamChildrenIndex)) = [];
        else
            stf(i).propVMAT.numOfBeamChildren = 0;
        end
        if isfield(stf(i).propVMAT,'beamSubChildrenGantryAngles')
            stf(i).propVMAT.beamSubChildrenGantryAngles(isnan(stf(i).propVMAT.beamSubChildrenGantryAngles)) = [];
            stf(i).propVMAT.beamSubChildrenIndex(isnan(stf(i).propVMAT.beamSubChildrenIndex)) = [];
        else
            stf(i).propVMAT.numOfBeamSubChildren = 0;
        end
    end
    
    if stf(i).propVMAT.DAOBeam
        if pln.propOpt.VMAToptions.continuousAperture
            
            stf(i).propVMAT.doseAngleDAO = ones(1,2);
            
            if sum(DAOAngleBorders == stf(i).propVMAT.DAOAngleBorders(2)) > 1
                %final dose angle is repeated
                %do not count twice in optimization
                stf(i).propVMAT.doseAngleDAO(2) = 0;
            end
        end
    else
        % dummy variable for doseAngleDAO
        stf(i).propVMAT.doseAngleDAO = [];
    end
    
    % determine the real last and next DAO indices
    % these are the indices whose DAOAngleBorders surround the
    % fluAngleBorders (1) and (2)
    lastDAOIndex_DAO    = find(abs(pln.propStf.fluGantryAngles(stf(i).propVMAT.lastDAOIndex) - pln.propStf.DAOGantryAngles) < 1e-8)+1;
    lastDAOIndex        = find(abs(pln.propStf.fluGantryAngles - pln.propStf.DAOGantryAngles(lastDAOIndex_DAO)) < 1e-8);
    while ~(stf(lastDAOIndex).propVMAT.DAOAngleBorders(1) <= stf(i).propVMAT.fluAngleBorders(1) && stf(i).propVMAT.fluAngleBorders(1) <= stf(lastDAOIndex).propVMAT.DAOAngleBorders(2)) && lastDAOIndex_DAO > 1
        
        lastDAOIndex_DAO    = lastDAOIndex_DAO-1;
        lastDAOIndex        = find(abs(pln.propStf.fluGantryAngles - pln.propStf.DAOGantryAngles(lastDAOIndex_DAO)) < 1e-8);
    end
    stf(i).propVMAT.lastDAOIndex = lastDAOIndex;
    
    nextDAOIndex_DAO    = find(abs(pln.propStf.fluGantryAngles(stf(i).propVMAT.nextDAOIndex) - pln.propStf.DAOGantryAngles) < 1e-8)-1;
    nextDAOIndex        = find(abs(pln.propStf.fluGantryAngles - pln.propStf.DAOGantryAngles(nextDAOIndex_DAO)) < 1e-8);
    while ~(stf(nextDAOIndex).propVMAT.DAOAngleBorders(1) <= stf(i).propVMAT.fluAngleBorders(2) && stf(i).propVMAT.fluAngleBorders(2) <= stf(nextDAOIndex).propVMAT.DAOAngleBorders(2)) && nextDAOIndex_DAO < numel(pln.propStf.DAOGantryAngles)
        
        nextDAOIndex_DAO    = nextDAOIndex_DAO+1;
        nextDAOIndex        = find(abs(pln.propStf.fluGantryAngles - pln.propStf.DAOGantryAngles(nextDAOIndex_DAO)) < 1e-8);
    end
    stf(i).propVMAT.nextDAOIndex = nextDAOIndex;
    
    % determine the real last and next dose indices
    % these are the indices whose doseAngleBorders surround the
    % fluAngleBorders (1) and (2)
    lastDoseIndex_dose  = find(abs(pln.propStf.fluGantryAngles(stf(i).propVMAT.lastDoseIndex) - pln.propStf.gantryAngles) < 1e-8)+1;
    lastDoseIndex       = find(abs(pln.propStf.fluGantryAngles - pln.propStf.gantryAngles(lastDoseIndex_dose)) < 1e-8);
    while ~(stf(lastDoseIndex).propVMAT.doseAngleBorders(1) <= stf(i).propVMAT.fluAngleBorders(1) && stf(i).propVMAT.fluAngleBorders(1) <= stf(lastDoseIndex).propVMAT.doseAngleBorders(2)) && lastDoseIndex_dose > 1
        
        lastDoseIndex_dose  = lastDoseIndex_dose-1;
        lastDoseIndex       = find(abs(pln.propStf.fluGantryAngles - pln.propStf.gantryAngles(lastDoseIndex_dose)) < 1e-8);
    end
    stf(i).propVMAT.lastDoseIndex = lastDoseIndex;
    
    nextDoseIndex_dose  = find(abs(pln.propStf.fluGantryAngles(stf(i).propVMAT.nextDoseIndex) - pln.propStf.gantryAngles) < 1e-8)-1;
    nextDoseIndex       = find(abs(pln.propStf.fluGantryAngles - pln.propStf.gantryAngles(nextDoseIndex_dose)) < 1e-8);
    while ~(stf(nextDoseIndex).propVMAT.doseAngleBorders(1) <= stf(i).propVMAT.fluAngleBorders(2) && stf(i).propVMAT.fluAngleBorders(2) <= stf(nextDoseIndex).propVMAT.doseAngleBorders(2)) && nextDoseIndex_dose < numel(pln.propStf.gantryAngles)
        
        nextDoseIndex_dose  = nextDoseIndex_dose+1;
        nextDoseIndex       = find(abs(pln.propStf.fluGantryAngles - pln.propStf.gantryAngles(nextDoseIndex_dose)) < 1e-8);
    end
    stf(i).propVMAT.nextDoseIndex = nextDoseIndex;
    
    %{
    % determine the current DAO index
    if stf(stf(i).propVMAT.lastDAOIndex).propVMAT.DAOAngleBorders(1) <= stf(i).gantryAngle && stf(i).gantryAngle <= stf(stf(i).propVMAT.lastDAOIndex).propVMAT.DAOAngleBorders(2)
        stf(i).propVMAT.currDAOIndex = stf(i).propVMAT.lastDAOIndex;
    elseif stf(stf(i).propVMAT.nextDAOIndex).propVMAT.DAOAngleBorders(1) <= stf(i).gantryAngle && stf(i).gantryAngle <= stf(stf(i).propVMAT.nextDAOIndex).propVMAT.DAOAngleBorders(2)
        stf(i).propVMAT.currDAOIndex = stf(i).propVMAT.nextDAOIndex;
    end
    
    % determine the current dose index
    if stf(stf(i).propVMAT.lastDoseIndex).propVMAT.doseAngleBorders(1) <= stf(i).gantryAngle && stf(i).gantryAngle <= stf(stf(i).propVMAT.lastDoseIndex).propVMAT.doseAngleBorders(2)
        stf(i).propVMAT.currDoseIndex = stf(i).propVMAT.lastDoseIndex;
    elseif stf(stf(i).propVMAT.nextDoseIndex).propVMAT.doseAngleBorders(1) <= stf(i).gantryAngle && stf(i).gantryAngle <= stf(stf(i).propVMAT.nextDoseIndex).propVMAT.doseAngleBorders(2)
        stf(i).propVMAT.currDoseIndex = stf(i).propVMAT.nextDoseIndex;
    end
    %}
    
    % calculate fractions for leaf position interpolation
    % formerly known as fracFromLastDAO_I, fracFromLastDAO_F
    stf(i).propVMAT.fracFromLastDAOI_leafI = (stf(stf(i).propVMAT.lastDAOIndex).propVMAT.DAOAngleBorders(2)-stf(i).propVMAT.fluAngleBorders(1))./stf(stf(i).propVMAT.lastDAOIndex).propVMAT.DAOAngleBordersDiff;
    stf(i).propVMAT.fracFromLastDAOF_leafI = 1-stf(i).propVMAT.fracFromLastDAOI_leafI;
    stf(i).propVMAT.fracFromNextDAOI_leafF = (stf(stf(i).propVMAT.nextDAOIndex).propVMAT.DAOAngleBorders(2)-stf(i).propVMAT.fluAngleBorders(2))./stf(stf(i).propVMAT.nextDAOIndex).propVMAT.DAOAngleBordersDiff;
    stf(i).propVMAT.fracFromNextDAOF_leafF = 1-stf(i).propVMAT.fracFromNextDAOI_leafF;
    
    % calculate fractions for dose calculation
    stf(i).propVMAT.fracToLastDose_arcI = (stf(stf(i).propVMAT.lastDoseIndex).propVMAT.doseAngleBorders(2)-stf(i).propVMAT.fluAngleBorders(1))./stf(i).propVMAT.fluAngleBorderCentreDiff(1);
    stf(i).propVMAT.fracToNextDose_arcF = (stf(i).propVMAT.fluAngleBorders(2)-stf(stf(i).propVMAT.nextDoseIndex).propVMAT.doseAngleBorders(1))./stf(i).propVMAT.fluAngleBorderCentreDiff(2);
    % fix instances where the fraction is outside the range [0 1]
    stf(i).propVMAT.fracToLastDose_arcI = round(max(min(stf(i).propVMAT.fracToLastDose_arcI,1),0),10);
    stf(i).propVMAT.fracToNextDose_arcF = round(max(min(stf(i).propVMAT.fracToNextDose_arcF,1),0),10);
    % the fraction to the next and last dose beams add up to 1
    stf(i).propVMAT.fracToNextDose_arcI = 1-stf(i).propVMAT.fracToLastDose_arcI;
    stf(i).propVMAT.fracToLastDose_arcF = 1-stf(i).propVMAT.fracToNextDose_arcF;
    
    
    % calculate MU rate fractions
    % technically only needed for non-DAO beams
    % formerly known as fracFromLastDAO
    stf(i).propVMAT.fracFromLastDAO_MU = (pln.propStf.fluGantryAngles(nextDAOIndex)-pln.propStf.fluGantryAngles(i))./(pln.propStf.fluGantryAngles(nextDAOIndex)-pln.propStf.fluGantryAngles(lastDAOIndex));
    % fix instances where this fraction is outside the range [0 1]
    % this can happen with fluence beams before the first / after the
    % last DAO beam
    stf(i).propVMAT.fracFromLastDAO_MU = max(min(stf(i).propVMAT.fracFromLastDAO_MU,1),0);
    % the fraction from the next and last DAO beams add up to 1
    stf(i).propVMAT.fracFromNextDAO_MU = 1-stf(i).propVMAT.fracFromLastDAO_MU;
    
    
    % calculate fractions for gantry rotation speed interpolation
    % technically only needed for non-DAO beams
    % formerly known as timeFracFromLastDAO
    stf(i).propVMAT.fracFromLastDAO_gantryRot = (stf(stf(i).propVMAT.lastDAOIndex).propVMAT.DAOAngleBorders(2)-stf(i).propVMAT.fluAngleBorders(1))./stf(i).propVMAT.fluAngleBordersDiff;
    % fix instances where the fraction is outside the range [0 1]
    % this can happen with fluence beams before the first / after the
    % last DAO beam
    stf(i).propVMAT.fracFromLastDAO_gantryRot = max(min(stf(i).propVMAT.fracFromLastDAO_gantryRot,1),0);
    % the fraction from the next and last DAO beams add up to 1
    stf(i).propVMAT.fracFromNextDAO_gantryRot = 1-stf(i).propVMAT.fracFromLastDAO_gantryRot;
    
    matRad_progress(i,length(stf));
end

end
