function jacob = matRad_daoJacobFunc(apertureInfoVec,dij,cst,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: jacobian function for direct aperture optimization
%
% call
%   jacob = matRad_daoJacobFunc(apertureInfoVec,apertureInfo,dij,cst,type)
%
% input
%   apertureInfoVec: aperture info vector
%   apertureInfo:    aperture info struct
%   dij:             dose influence matrix
%   cst:             matRad cst struct
%   options:         option struct defining the type of optimization
%
% output
%   jacob:           jacobian of constraint function
%
% References
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

% jacobian of the dao constraints

% row indices
i = repmat(1:apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases,1,2);
% column indices
j = [(apertureInfo.totalNumOfShapes*apertureInfo.numPhases+1):(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs)*apertureInfo.numPhases ...
    ((apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs)*apertureInfo.numPhases+1):(apertureInfo.totalNumOfShapes+2*apertureInfo.totalNumOfLeafPairs)*apertureInfo.numPhases];

% -1 for left leaves, 1 for right leaves
s = [-1*ones(1,apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases) ones(1,apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases)];

jacob_dao = sparse(i,j,s, ...
    apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases, ...
    numel(apertureInfoVec), ...
    2*apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases);

% compute jacobian of dosimetric constrainst

% dosimetric jacobian in bixel space
jacob_dos_bixel = matRad_jacobFuncWrapper(apertureInfo.bixelWeights,dij,cst,options);

% is this necessary?
if ~isempty(jacob_dos_bixel{1})
    
    numOfConstraints = size(jacob_dos_bixel{1},1);
    
    i_sparse = 1:numOfConstraints;
    i_sparse = kron(i_sparse,ones(1,numel(apertureInfoVec)));
    
    j_sparse = 1:numel(apertureInfoVec);
    j_sparse = repmat(j_sparse,1,numOfConstraints);
    
    jacobSparseVec = zeros(numOfConstraints*size(apertureInfoVec,1),1);
    
    if apertureInfo.runVMAT
        
        jacob_dos = sparse(numOfConstraints,numel(apertureInfoVec));
        
        for phase = 1:apertureInfo.numPhases
            
            % use the Jacobian calculated in daoVec2ApertureInfo.
            % should also do this for non-VMAT
            jacob_dos = jacob_dos+jacob_dos_bixel{phase}*apertureInfo.bixelJApVec{phase}';
        end
        
    else
        
        % 1. calculate jacobian for aperture weights
        % loop over all beams
        conOffset = 0;
        for i = 1:numel(apertureInfo.beam)
            
            % get used bixels in beam
            ix = ~isnan(apertureInfo.beam(i).bixelIndMap);
            
            % loop over all shapes and add up the gradients x openingFrac for this shape
            for j = 1:apertureInfo.beam(i).numOfShapes
                
                jacobSparseVec(conOffset+j == j_sparse) = jacob_dos_bixel{1}(:,apertureInfo.beam(i).bixelIndMap(ix)) ...
                    * apertureInfo.beam(i).shape(j).shapeMap(ix)./apertureInfo.beam(i).shape(j).jacobiScale;
            end
            
            % increment offset
            conOffset = conOffset + apertureInfo.beam(i).numOfShapes;
            
        end
        
        % 2. find corresponding bixel to the leaf Positions and aperture
        % weights to calculate the jacobian
        indInSparseVec = repmat(apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2,1,numOfConstraints) ...
            +kron((0:numOfConstraints-1)*numel(apertureInfoVec),ones(1,apertureInfo.totalNumOfLeafPairs*2));
        
        jacobSparseVec(indInSparseVec) = ...
            reshape(transpose(( ones(numOfConstraints,1) * apertureInfoVec(apertureInfo.mappingMx(apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2,2))' ) ...
            .* jacob_dos_bixel{1}(:,apertureInfo.bixelIndices(apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)) ./ ...
            (ones(numOfConstraints,1) * (apertureInfo.bixelWidth.*apertureInfo.jacobiScale(apertureInfo.mappingMx(apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2,2)))')),[],1);
        
        
        % correct the sign for the left leaf positions
        indInSparseVec = repmat(apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs,1,numOfConstraints) ...
            +kron((0:numOfConstraints-1)*numel(apertureInfoVec),ones(1,apertureInfo.totalNumOfLeafPairs));
        
        jacobSparseVec(indInSparseVec) = -jacobSparseVec(indInSparseVec);
        
        jacob_dos = sparse(i_sparse,j_sparse,jacobSparseVec,numOfConstraints,numel(apertureInfoVec));
    end
    
else
    jacob_dos = sparse(0,0);
end

if ~apertureInfo.runVMAT
    % concatenate
    jacob = [jacob_dao; jacob_dos];
else
    
    % get timeFacCurr, _rep
    timeFacCurr     = [apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).timeFacCurr]';
    timeFacCurr_rep = repmat(timeFacCurr,apertureInfo.numPhases,1);
    
    if apertureInfo.propVMAT.continuousAperture
        
        % set up
        n = apertureInfo.beam(1).numOfActiveLeafPairs;
        indInSparseVec  = (1:n);
        indInConVec     = (1:n);
        
        % sparse matrix
        if apertureInfo.propVMAT.fixedGantrySpeed
            numElem = n.*apertureInfo.propVMAT.numLeafSpeedConstraint.*8;
        else
            numElem = n.*(apertureInfo.propVMAT.numLeafSpeedConstraintDAO.*10 + (apertureInfo.propVMAT.numLeafSpeedConstraint-apertureInfo.propVMAT.numLeafSpeedConstraintDAO).*12);
        end
        i_sparse    = zeros(numElem,1);
        j_sparse    = zeros(numElem,1);
        s_sparse    = zeros(numElem,1);
        
        for i = 1:numel(apertureInfo.beam)
            % loop over beams
            
            % extract time spent in fluence angle arc
            tFluBorderAngle = apertureInfo.beam(i).time;
            
            % extract factors
            fracFromLastDAOI_leafI = apertureInfo.propVMAT.beam(i).fracFromLastDAOI_leafI;
            fracFromLastDAOF_leafI = apertureInfo.propVMAT.beam(i).fracFromLastDAOF_leafI;
            fracFromNextDAOI_leafF = apertureInfo.propVMAT.beam(i).fracFromNextDAOI_leafF;
            fracFromNextDAOF_leafF = apertureInfo.propVMAT.beam(i).fracFromNextDAOF_leafF;
            
            if apertureInfo.propVMAT.beam(i).DAOBeam
                jacobT_curr = apertureInfo.propVMAT.jacobT(apertureInfo.propVMAT.beam(i).DAOIndex,i);
            else
                jacobT_last = apertureInfo.propVMAT.jacobT(apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).DAOIndex,i);
                jacobT_next = apertureInfo.propVMAT.jacobT(apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).DAOIndex,i);
            end
            
            for phase_I = 1:apertureInfo.numPhases
                % loop over initial phases
                
                % find allowable transitions
                transMask                   = apertureInfo.propVMAT.beam(i).transMask(phase_I,:);
                transMask(transMask == 0)   = [];
                
                for phase_F = transMask
                    % loop over possible final phases
                    
                    % get vector indices
                    vectorIx_L_lastDAOI = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_I}.vectorOffset(1) + ((1:n)-1);
                    vectorIx_L_lastDAOF = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase_I}.vectorOffset(2) + ((1:n)-1);
                    vectorIx_L_nextDAOI = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_F}.vectorOffset(1) + ((1:n)-1);
                    vectorIx_L_nextDAOF = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase_F}.vectorOffset(2) + ((1:n)-1);
                    
                    vectorIx_R_lastDAOI = vectorIx_L_lastDAOI+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
                    vectorIx_R_lastDAOF = vectorIx_L_lastDAOF+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
                    vectorIx_R_nextDAOI = vectorIx_L_nextDAOI+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
                    vectorIx_R_nextDAOF = vectorIx_L_nextDAOF+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
                    
                    % extract leaf positions
                    leftLeafPos_I   = apertureInfo.beam(i).shape{phase_I}.leftLeafPos_I;
                    leftLeafPos_F   = apertureInfo.beam(i).shape{phase_F}.leftLeafPos_F;
                    rightLeafPos_I  = apertureInfo.beam(i).shape{phase_I}.rightLeafPos_I;
                    rightLeafPos_F  = apertureInfo.beam(i).shape{phase_F}.rightLeafPos_F;
                    
                    % calc diffs
                    leftLeafDiff    = leftLeafPos_F-leftLeafPos_I;
                    rightLeafDiff   = rightLeafPos_F-rightLeafPos_I;
                    
                    % calc jacobs
                    
                    % wrt last initial left (optimization vector)
                    i_sparse(indInSparseVec)    = indInConVec;
                    j_sparse(indInSparseVec)    = vectorIx_L_lastDAOI;
                    s_sparse(indInSparseVec)    = -fracFromLastDAOI_leafI.*sign(leftLeafDiff)./tFluBorderAngle;
                    indInSparseVec              = indInSparseVec+n;
                    
                    % wrt last final left (optimization vector)
                    i_sparse(indInSparseVec)    = indInConVec;
                    j_sparse(indInSparseVec)    = vectorIx_L_lastDAOF;
                    s_sparse(indInSparseVec)    = -fracFromLastDAOF_leafI.*sign(leftLeafDiff)./tFluBorderAngle;
                    indInSparseVec              = indInSparseVec+n;
                    
                    % wrt next initial left (optimization vector)
                    i_sparse(indInSparseVec)    = indInConVec;
                    j_sparse(indInSparseVec)    = vectorIx_L_nextDAOI;
                    s_sparse(indInSparseVec)    = fracFromNextDAOI_leafF.*sign(leftLeafDiff)./tFluBorderAngle;
                    indInSparseVec              = indInSparseVec+n;
                    
                    % wrt next final left (optimization vector)
                    i_sparse(indInSparseVec)    = indInConVec;
                    j_sparse(indInSparseVec)    = vectorIx_L_nextDAOF;
                    s_sparse(indInSparseVec)    = fracFromNextDAOF_leafF.*sign(leftLeafDiff)./tFluBorderAngle;
                    indInSparseVec              = indInSparseVec+n;
                    
                    % wrt last initial right (optimization vector)
                    i_sparse(indInSparseVec)    = indInConVec+apertureInfo.propVMAT.numLeafSpeedConstraint*apertureInfo.beam(1).numOfActiveLeafPairs;
                    j_sparse(indInSparseVec)    = vectorIx_R_lastDAOI;
                    s_sparse(indInSparseVec)    = -fracFromLastDAOI_leafI.*sign(rightLeafDiff)./tFluBorderAngle;
                    indInSparseVec              = indInSparseVec+n;
                    
                    % wrt last final right (optimization vector)
                    i_sparse(indInSparseVec)    = indInConVec+apertureInfo.propVMAT.numLeafSpeedConstraint*apertureInfo.beam(1).numOfActiveLeafPairs;
                    j_sparse(indInSparseVec)    = vectorIx_R_lastDAOF;
                    s_sparse(indInSparseVec)    = -fracFromLastDAOF_leafI.*sign(rightLeafDiff)./tFluBorderAngle;
                    indInSparseVec              = indInSparseVec+n;
                    
                    % wrt next initial right (optimization vector)
                    i_sparse(indInSparseVec)    = indInConVec+apertureInfo.propVMAT.numLeafSpeedConstraint*apertureInfo.beam(1).numOfActiveLeafPairs;
                    j_sparse(indInSparseVec)    = vectorIx_R_nextDAOI;
                    s_sparse(indInSparseVec)    = fracFromNextDAOI_leafF.*sign(rightLeafDiff)./tFluBorderAngle;
                    indInSparseVec              = indInSparseVec+n;
                    
                    % wrt next final right (optimization vector)
                    i_sparse(indInSparseVec)    = indInConVec+apertureInfo.propVMAT.numLeafSpeedConstraint*apertureInfo.beam(1).numOfActiveLeafPairs;
                    j_sparse(indInSparseVec)    = vectorIx_R_nextDAOF;
                    s_sparse(indInSparseVec)    = fracFromNextDAOF_leafF.*sign(rightLeafDiff)./tFluBorderAngle;
                    indInSparseVec              = indInSparseVec+n;
                    
                    if ~apertureInfo.propVMAT.fixedGantrySpeed
                        % only do time Jacobian if we are doing variable
                        % gantry speed
                        if apertureInfo.propVMAT.beam(i).DAOBeam
                            % wrt time (left, then right)
                            % it's a DAO beam, so speed only depends on its own
                            % time
                            i_sparse(indInSparseVec)    = indInConVec;
                            j_sparse(indInSparseVec)    = apertureInfo.propVMAT.beam(i).timeInd;
                            s_sparse(indInSparseVec)    = -abs(leftLeafDiff).*jacobT_curr./(tFluBorderAngle.^2);
                            indInSparseVec              = indInSparseVec+n;
                            
                            i_sparse(indInSparseVec)    = indInConVec+apertureInfo.propVMAT.numLeafSpeedConstraint*apertureInfo.beam(1).numOfActiveLeafPairs;
                            j_sparse(indInSparseVec)    = apertureInfo.propVMAT.beam(i).timeInd;
                            s_sparse(indInSparseVec)    = -abs(rightLeafDiff).*jacobT_curr./(tFluBorderAngle.^2);
                            indInSparseVec              = indInSparseVec+n;
                        else
                            % this is not a DAO beam, so the speed will
                            % depend on the last and next DAO indices
                            
                            % first do last (left, then right)
                            i_sparse(indInSparseVec)    = indInConVec;
                            j_sparse(indInSparseVec)    = apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).timeInd;
                            s_sparse(indInSparseVec)    = -abs(leftLeafDiff).*jacobT_last./(tFluBorderAngle.^2);
                            indInSparseVec              = indInSparseVec+n;
                            
                            i_sparse(indInSparseVec)    = indInConVec+apertureInfo.propVMAT.numLeafSpeedConstraint*apertureInfo.beam(1).numOfActiveLeafPairs;
                            j_sparse(indInSparseVec)    = apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).timeInd;
                            s_sparse(indInSparseVec)    = -abs(rightLeafDiff).*jacobT_last./(tFluBorderAngle.^2);
                            indInSparseVec              = indInSparseVec+n;
                            
                            if apertureInfo.propVMAT.beam(i).nextDAOIndex ~= apertureInfo.propVMAT.beam(i).lastDAOIndex
                                % now do next only if different from last 
                                % (left, then right)
                                i_sparse(indInSparseVec)    = indInConVec;
                                j_sparse(indInSparseVec)    = apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).timeInd;
                                s_sparse(indInSparseVec)    = -abs(leftLeafDiff).*jacobT_next./(tFluBorderAngle.^2);
                                indInSparseVec              = indInSparseVec+n;
                                
                                i_sparse(indInSparseVec)    = indInConVec+apertureInfo.propVMAT.numLeafSpeedConstraint*apertureInfo.beam(1).numOfActiveLeafPairs;
                                j_sparse(indInSparseVec)    = apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).timeInd;
                                s_sparse(indInSparseVec)    = -abs(rightLeafDiff).*jacobT_next./(tFluBorderAngle.^2);
                                indInSparseVec              = indInSparseVec+n;
                            end
                        end
                    end
                    
                    % update offset
                    indInConVec = indInConVec+n;
                end
            end
        end
        
        i_sparse(indInSparseVec(1):end) = [];
        j_sparse(indInSparseVec(1):end) = [];
        s_sparse(indInSparseVec(1):end) = [];
        
        jacob_lfspd = sparse(i_sparse,j_sparse,s_sparse,2*apertureInfo.beam(1).numOfActiveLeafPairs*apertureInfo.propVMAT.numLeafSpeedConstraint,numel(apertureInfoVec));
        
    else
        
        % values of times spent in an arc surrounding the optimized angles
        % (full arc)
        timeFluBorderAngles     = [apertureInfo.beam([apertureInfo.propVMAT.beam.DAOBeam]).time]';
        timeDAOBorderAngles     = timeFluBorderAngles./timeFacCurr;
        
        % get index values for the jacobian
        % variable index
        % value of constraints for leaves
        leftLeafPos  = apertureInfoVec((1:apertureInfo.totalNumOfLeafPairs)+apertureInfo.totalNumOfShapes);
        rightLeafPos = apertureInfoVec(1+apertureInfo.totalNumOfLeafPairs+apertureInfo.totalNumOfShapes:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2);
        
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
        
        currentLeftLeafInd = (apertureInfo.totalNumOfShapes+1):(apertureInfo.totalNumOfShapes+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeBNOptAngles));
        currentRightLeafInd = (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+1):(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeBNOptAngles));
        nextLeftLeafInd = (apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+1):(apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeBNOptAngles));
        nextRightLeafInd = (apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+1):(apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeBNOptAngles));
        leftTimeInd = kron(j,ones(1,apertureInfo.beam(1).numOfActiveLeafPairs))+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2;
        rightTimeInd = kron(j,ones(1,apertureInfo.beam(1).numOfActiveLeafPairs))+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2;
        %leftTimeInd = repelem(j,apertureInfo.beam(1).numOfActiveLeafPairs)+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2;
        %rightTimeInd = repelem(j,apertureInfo.beam(1).numOfActiveLeafPairs)+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2;
        % constraint index
        constraintInd = 1:2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeBNOptAngles);
        
        
        % jacobian of the leafspeed constraint
        i = repmat((i'-1)*apertureInfo.beam(1).numOfActiveLeafPairs,1,apertureInfo.beam(1).numOfActiveLeafPairs)+repmat(1:apertureInfo.beam(1).numOfActiveLeafPairs,2*numel(timeBNOptAngles),1);
        i = reshape([i' i'+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeBNOptAngles)],1,[]);
        
        i = [repmat(constraintInd,1,2) i];
        j = [currentLeftLeafInd currentRightLeafInd nextLeftLeafInd nextRightLeafInd leftTimeInd rightTimeInd];
        % first do jacob wrt current leaf position (left, right), then next leaf
        % position (left, right), then time (left, right)
        j_lfspd_cur = -reshape([sign(diff(reshape(leftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2)) ...
            sign(diff(reshape(rightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2))]./ ...
            repmat(timeBNOptAngles',apertureInfo.beam(1).numOfActiveLeafPairs,2),2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeBNOptAngles),1);
        
        j_lfspd_nxt = reshape([sign(diff(reshape(leftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2)) ...
            sign(diff(reshape(rightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2))]./ ...
            repmat(timeBNOptAngles',apertureInfo.beam(1).numOfActiveLeafPairs,2),2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeBNOptAngles),1);
        
        j_lfspd_t = -reshape([kron(abs(diff(reshape(leftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2)),ones(1,2)).*repmat(timeFac',apertureInfo.beam(1).numOfActiveLeafPairs,1) ...
            kron(abs(diff(reshape(rightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2)),ones(1,2)).*repmat(timeFac',apertureInfo.beam(1).numOfActiveLeafPairs,1)]./ ...
            repmat(kron((timeBNOptAngles.^2)',ones(1,2)),apertureInfo.beam(1).numOfActiveLeafPairs,2),[],1);
        
        s = [j_lfspd_cur; j_lfspd_nxt; j_lfspd_t];
        
        jacob_lfspd = sparse(i,j,s,2*apertureInfo.beam(1).numOfActiveLeafPairs*(apertureInfo.totalNumOfShapes-1),numel(apertureInfoVec));
    end
    
    % jacobian of the doserate constraint
    % values of doserate (MU/sec) between optimized gantry angles
    weights = apertureInfoVec(1:(apertureInfo.totalNumOfShapes*apertureInfo.numPhases))./apertureInfo.jacobiScale;
    
    % values of times spent in an arc surrounding the optimized angles
    % (dose influence arc)
    timeFluBorderAngles     = [apertureInfo.beam([apertureInfo.propVMAT.beam.DAOBeam]).time]';
    timeFluBorderAngles_rep = repmat(timeFluBorderAngles,apertureInfo.numPhases,1);
    
    % first do jacob wrt weights, then wrt times (if applicable)
    if apertureInfo.propVMAT.fixedGantrySpeed
        i = 1:(apertureInfo.totalNumOfShapes*apertureInfo.numPhases);
        j = 1:(apertureInfo.totalNumOfShapes*apertureInfo.numPhases);
        s = apertureInfo.weightToMU./(timeFluBorderAngles_rep.*apertureInfo.jacobiScale);
    else
        i = repmat(1:(apertureInfo.totalNumOfShapes*apertureInfo.numPhases),1,2);
        j = [1:(apertureInfo.totalNumOfShapes*apertureInfo.numPhases) ...
            repmat(((apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases+1):((apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases+apertureInfo.totalNumOfShapes),1,apertureInfo.numPhases)];
        s = [apertureInfo.weightToMU./(timeFluBorderAngles_rep.*apertureInfo.jacobiScale); -apertureInfo.weightToMU.*weights.*timeFacCurr_rep./(timeFluBorderAngles_rep.^2)];
    end
    
    jacob_dosrt = sparse(i,j,s,apertureInfo.totalNumOfShapes*apertureInfo.numPhases,numel(apertureInfoVec));
    
    % concatenate
    jacob = [jacob_dao; jacob_lfspd; jacob_dosrt; jacob_dos];
end


