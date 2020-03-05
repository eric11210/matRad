function apertureInfo = matRad_maxLeafSpeed(apertureInfo)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad calculation of maximum leaf speed
%
% call
%   apertureInfo = matRad_maxLeafSpeed(apertureInfo)
%
% input
%   apertureInfo:   aperture info struct
%
% output
%   apertureInfo:   aperture info struct
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


apertureInfoVec = apertureInfo.apertureVector;

% find values of leaf speeds of optimized gantry angles
if apertureInfo.propVMAT.continuousAperture
    % Using the dynamic fluence calculation, we have the leaf positions in
    % the vector be the leaf positions at the borders of the DAO arcs
    % (for optimized angles only).
    % Therefore we must also use the times between the borders of the Dij
    % arc (for optimized angles only).
    
    % prep
    leftLeafDiff    = zeros(apertureInfo.propVMAT.numLeafSpeedConstraint*apertureInfo.beam(1).numOfActiveLeafPairs,1);
    rightLeafDiff   = zeros(apertureInfo.propVMAT.numLeafSpeedConstraint*apertureInfo.beam(1).numOfActiveLeafPairs,1);
    tVec            = zeros(apertureInfo.propVMAT.numLeafSpeedConstraint*apertureInfo.beam(1).numOfActiveLeafPairs,1);
    maxLeafSpeed    = zeros(1,numel(apertureInfo.beam));
    
    offset      = 0;
    
    for i = 1:numel(apertureInfo.beam)
        % loop over beams
        n = apertureInfo.beam(i).numOfActiveLeafPairs;
        
        % extract time spent in fluence angle arc
        tFluBorderAngle = apertureInfo.beam(i).time;
        
        for phase_I = 1:apertureInfo.numPhases
            % loop over initial phases
            
            % find allowable transitions
            transMask                   = apertureInfo.propVMAT.beam(i).transMask(phase_I,:);
            transMask(transMask == 0)   = [];
            
            for phase_F = transMask
                % loop over possible final phases
                
                % extract leaf positions
                leftLeafPos_I   = apertureInfo.beam(i).shape{phase_I}.leftLeafPos_I;
                leftLeafPos_F   = apertureInfo.beam(i).shape{phase_F}.leftLeafPos_F;
                rightLeafPos_I  = apertureInfo.beam(i).shape{phase_I}.rightLeafPos_I;
                rightLeafPos_F  = apertureInfo.beam(i).shape{phase_F}.rightLeafPos_F;
                
                % determine indices
                indInDiffVec = offset+(1:n);
                
                % insert differences, time
                leftLeafDiff(indInDiffVec)  = abs(leftLeafPos_F-leftLeafPos_I);
                rightLeafDiff(indInDiffVec) = abs(rightLeafPos_F-rightLeafPos_I);
                tVec(indInDiffVec)          = tFluBorderAngle;
                
                % get max speed
                leftLeafSpeed = abs(leftLeafPos_F-leftLeafPos_I)./tFluBorderAngle;
                rightLeafSpeed = abs(leftLeafPos_F-leftLeafPos_I)./tFluBorderAngle;
                maxLeafSpeed_temp = max([leftLeafSpeed; rightLeafSpeed]);
                
                % update max speed
                if maxLeafSpeed_temp > maxLeafSpeed(i)
                    maxLeafSpeed(i) = maxLeafSpeed_temp;
                end
                
                % update offset
                offset = offset+n;
            end
        end
    end
else
    % values of time differences of optimized gantry angles
    timeDAOBorderAngles = [apertureInfo.beam([apertureInfo.propVMAT.beam.DAOBeam]).time]./[apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).timeFacCurr];
    
    % value of constraints for leaves
    %leftLeafPos  = apertureInfoVec([1:apertureInfo.totalNumOfLeafPairs]+apertureInfo.totalNumOfShapes);
    %rightLeafPos = apertureInfoVec(1+apertureInfo.totalNumOfLeafPairs+apertureInfo.totalNumOfShapes:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2);
    leftLeafPos  = apertureInfoVec((1:(apertureInfo.numPhases*apertureInfo.totalNumOfLeafPairs))+apertureInfo.numPhases*apertureInfo.totalNumOfShapes);
    rightLeafPos = apertureInfoVec(1+(apertureInfo.totalNumOfLeafPairs+apertureInfo.totalNumOfShapes)*apertureInfo.numPhases:(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases);
    
    % Using the static fluence calculation, we have the leaf positions in
    % the vector be the leaf positions at the centre of the Dij arcs (for optimized angles only).
    % Therefore we must use the times between the centres of the Dij arcs (for optimized angles only).
    i = sort(repmat(1:(apertureInfo.totalNumOfShapes-1),1,2));
    j = sort(repmat(1:apertureInfo.totalNumOfShapes,1,2));
    j(1) = [];
    j(end) = [];
    
    timeFac = [apertureInfo.propVMAT.beam.timeFac]';
    timeFac(1) = [];
    timeFac(end) = [];
    %timeFac(timeFac == 0) = [];
    
    timeFacMatrix = sparse(i,j,timeFac,(apertureInfo.totalNumOfShapes-1),apertureInfo.totalNumOfShapes);
    timeBNOptAngles = timeFacMatrix*timeDAOBorderAngles;
    
    leftLeafSpeed = abs(diff(reshape(leftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,[]),1,2))./repmat(timeBNOptAngles',apertureInfo.beam(1).numOfActiveLeafPairs,1);
    rightLeafSpeed = abs(diff(reshape(rightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,[]),1,2))./repmat(timeBNOptAngles',apertureInfo.beam(1).numOfActiveLeafPairs,1);
    
    % values of max leaf speeds
    leftMaxLeafSpeed = max(leftLeafSpeed,[],1);
    rightMaxLeafSpeed = max(rightLeafSpeed,[],1);
    maxLeafSpeed = max([leftMaxLeafSpeed; rightMaxLeafSpeed],[],1);
end


% enter into apertureInfo
maxMaxLeafSpeed = 0;
for i = 1:size(apertureInfo.beam,2)
    if apertureInfo.propVMAT.continuousAperture
        % for dynamic, we take the max leaf speed to be the actual leaf
        % speed
        
        apertureInfo.beam(i).maxLeafSpeed = maxLeafSpeed(i);
        if apertureInfo.beam(i).maxLeafSpeed >= maxMaxLeafSpeed
            maxMaxLeafSpeed = apertureInfo.beam(i).maxLeafSpeed;
        end
    else
        % for static, we take the max leaf speed to be the max leaf
        % of two speeds, one being the speed in the first half-arc, the
        % second being the speed in the second half-arc (these will be
        % different in general)
        
        if i == 1
            apertureInfo.beam(i).maxLeafSpeed = maxLeafSpeed(i);
        elseif i == apertureInfo.totalNumOfShapes
            apertureInfo.beam(i).maxLeafSpeed = maxLeafSpeed(i-1);
        else
            %apertureInfo.beam(i).maxLeafSpeed = maxLeafSpeed(l-1)*apertureInfo.beam(i).timeFac(1)+maxLeafSpeed(l)*apertureInfo.beam(i).timeFac(2);
            apertureInfo.beam(i).maxLeafSpeed = max(maxLeafSpeed(i-1),maxLeafSpeed(i));
        end
        
        
        if i < apertureInfo.totalNumOfShapes && maxLeafSpeed(i) >= maxMaxLeafSpeed
            maxMaxLeafSpeed = maxLeafSpeed(i);
        end
    end
end

apertureInfo.maxLeafSpeed = maxMaxLeafSpeed;

