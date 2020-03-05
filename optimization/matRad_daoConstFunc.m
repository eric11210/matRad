function c = matRad_daoConstFunc(apertureInfoVec,dij,cst,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: constraint function for direct aperture optimization
%
% call
%   c = matRad_daoObjFunc(apertueInfoVec,dij,cst)
%
% input
%   apertueInfoVec: aperture info vector
%   apertureInfo:   aperture info struct
%   dij:            dose influence matrix
%   cst:            matRad cst struct
%   options:        option struct defining the type of optimization
%
% output
%   c:              value of constraints
%
% Reference
%   [1] http://www.sciencedirect.com/science/article/pii/S0958394701000577
%   [2] http://www.sciencedirect.com/science/article/pii/S0360301601025858
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

% read in the global apertureInfo and apertureVector variables
global matRad_global_apertureInfo;
% update apertureInfo from the global variable
apertureInfo = matRad_global_apertureInfo;

% also read in the global recalc variance variable
global matRad_global_recalcVar;

% update apertureInfo, bixel weight vector an mapping of leafes to bixels
if ~isequal(apertureInfoVec,apertureInfo.apertureVector)
    apertureInfo = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfoVec);
    matRad_global_apertureInfo = apertureInfo;
    
    % recalculate the variance
    matRad_global_recalcVar = true;
end

% value of constraints for leaves
leftLeafPos  = apertureInfoVec((1:apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases)+apertureInfo.totalNumOfShapes*apertureInfo.numPhases);
rightLeafPos = apertureInfoVec((1+(apertureInfo.totalNumOfLeafPairs+apertureInfo.totalNumOfShapes)*apertureInfo.numPhases):(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases);
c_dao        = rightLeafPos - leftLeafPos;

% bixel based objective function calculation
c_dos = matRad_constFuncWrapper(apertureInfo.bixelWeights,dij,cst,options);

if ~apertureInfo.runVMAT
    
    % concatenate
    c = [c_dao; c_dos];
else
    
    if apertureInfo.propVMAT.continuousAperture
        % Using the dynamic fluence calculation, we have the leaf positions in
        % the vector be the leaf positions at the borders of the DAO arcs
        % (for optimized angles only).
        % Therefore we must also use the times between the borders of the Dij
        % arc (for optimized angles only).
        
        % prep
        leftLeafSpeed   = zeros(apertureInfo.propVMAT.numLeafSpeedConstraint*apertureInfo.beam(1).numOfActiveLeafPairs,1);
        rightLeafSpeed  = zeros(apertureInfo.propVMAT.numLeafSpeedConstraint*apertureInfo.beam(1).numOfActiveLeafPairs,1);
        
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
                    indInConVec = offset+(1:n);
                    
                    % calc speeds
                    leftLeafSpeed(indInConVec)      = abs(leftLeafPos_F-leftLeafPos_I)./tFluBorderAngle;
                    rightLeafSpeed(indInConVec)     = abs(rightLeafPos_F-rightLeafPos_I)./tFluBorderAngle;
                    
                    % update offset
                    offset = offset+n;
                end
            end
        end
        
        c_lfspd = [leftLeafSpeed; rightLeafSpeed];
    else
        
        % values of times spent in an arc surrounding the optimized angles
        % (full arc)
        timeFluBorderAngles     = [apertureInfo.beam([apertureInfo.propVMAT.beam.DAOBeam]).time]';
        timeFacCurr             = [apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).timeFacCurr]';
        timeDAOBorderAngles     = timeFluBorderAngles./timeFacCurr;
        
        i = sort(repmat(1:(apertureInfo.totalNumOfShapes-1),1,2));
        j = sort(repmat(1:apertureInfo.totalNumOfShapes,1,2));
        j(1) = [];
        j(end) = [];
        
        timeFac = [apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).timeFac]';
        timeFac(1) = [];
        timeFac(end) = [];
        %timeFac(timeFac == 0) = [];
        
        timeFacMatrix = sparse(i,j,timeFac,(apertureInfo.totalNumOfShapes-1),apertureInfo.totalNumOfShapes);
        timeBNOptAngles = timeFacMatrix*timeDAOBorderAngles;
        
        % values of average leaf speeds of optimized gantry angles
        c_lfspd = reshape([abs(diff(reshape(leftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2)) ...
            abs(diff(reshape(rightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2))]./ ...
            repmat(timeBNOptAngles',apertureInfo.beam(1).numOfActiveLeafPairs,2),2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeBNOptAngles),1);
    end
    
    % values of times spent in an arc surrounding the optimized angles
    % (dose influence arc)
    timeFluBorderAngles        = [apertureInfo.beam([apertureInfo.propVMAT.beam.DAOBeam]).time]';
    timeDoseBorderAngles_rep    = repmat(timeFluBorderAngles,apertureInfo.numPhases,1);
    
    % values of doserate (MU/sec) in an arc surrounding the optimized angles
    weights = apertureInfoVec(1:(apertureInfo.totalNumOfShapes*apertureInfo.numPhases))./apertureInfo.jacobiScale;
    c_dosrt = apertureInfo.weightToMU.*weights./timeDoseBorderAngles_rep;
    
    % concatenate
    c = [c_dao; c_lfspd; c_dosrt; c_dos];
    
end

