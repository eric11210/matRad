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

% update apertureInfo, bixel weight vector an mapping of leafes to bixels
if ~isequal(apertureInfoVec,apertureInfo.apertureVector)
    apertureInfo = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfoVec);
    matRad_global_apertureInfo = apertureInfo;
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
        offset = 0;
        for i = 1:numel(apertureInfo.beam);
            
            % get used bixels in beam
            ix = ~isnan(apertureInfo.beam(i).bixelIndMap);
            
            % loop over all shapes and add up the gradients x openingFrac for this shape
            for j = 1:apertureInfo.beam(i).numOfShapes
                
                jacobSparseVec(offset+j == j_sparse) = jacob_dos_bixel{1}(:,apertureInfo.beam(i).bixelIndMap(ix)) ...
                    * apertureInfo.beam(i).shape(j).shapeMap(ix)./apertureInfo.beam(i).shape(j).jacobiScale;
            end
            
            % increment offset
            offset = offset + apertureInfo.beam(i).numOfShapes;
            
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
    
    if apertureInfo.propVMAT.continuousAperture
        
        % value of constraints for leaves
        leftLeafPos  = apertureInfoVec([1:apertureInfo.totalNumOfLeafPairs]+apertureInfo.totalNumOfShapes);
        rightLeafPos = apertureInfoVec(1+apertureInfo.totalNumOfLeafPairs+apertureInfo.totalNumOfShapes:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2);
        
        % values of times spent in an arc surrounding the optimized angles (full
        % arc/dose influence arc)
        timeOptBorderAngles = apertureInfoVec((1+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2):(apertureInfo.totalNumOfShapes*2+apertureInfo.totalNumOfLeafPairs*2));
        timeFacCurr = [apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).timeFacCurr]';
        timeDoseBorderAngles = timeOptBorderAngles.*timeFacCurr;
        
        %indices for the leaf speed Jacobian
        initialLeftLeafInd = [apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).initialLeftLeafInd];
        initialRightLeafInd = [apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).initialRightLeafInd];
        finalLeftLeafInd = [apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).finalLeftLeafInd];
        finalRightLeafInd = [apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).finalRightLeafInd];
        leftTimeInd = [apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).leftTimeInd];
        rightTimeInd = [apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).rightTimeInd];
        % constraint index
        constraintInd = 1:2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeDoseBorderAngles);
        
        
        % jacobian of the leafspeed constraint
        i = repmat(constraintInd,1,3);
        j = [initialLeftLeafInd initialRightLeafInd finalLeftLeafInd finalRightLeafInd leftTimeInd rightTimeInd];
        % first do jacob wrt initial leaf position (left, right), then final leaf
        % position (left, right), then time (left, right)
        
        %%%% first do jacob wrt prev leaf positions, then next
        
        leftLeafDiff = diff(reshape(leftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,[]),1,2);
        rightLeafDiff = diff(reshape(rightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,[]),1,2);
        
        leftLeafDiff = reshape(leftLeafDiff(repmat([apertureInfo.propVMAT.beam.DAOBeam],apertureInfo.beam(1).numOfActiveLeafPairs,1)),apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes);
        rightLeafDiff = reshape(rightLeafDiff(repmat([apertureInfo.propVMAT.beam.DAOBeam],apertureInfo.beam(1).numOfActiveLeafPairs,1)),apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes);
        
        j_lfspd_cur = -reshape(sign([leftLeafDiff rightLeafDiff])./ ...
            repmat(timeDoseBorderAngles',apertureInfo.beam(1).numOfActiveLeafPairs,2),2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeDoseBorderAngles),1);
        
        j_lfspd_nxt = reshape(sign([leftLeafDiff rightLeafDiff])./ ...
            repmat(timeDoseBorderAngles',apertureInfo.beam(1).numOfActiveLeafPairs,2),2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeDoseBorderAngles),1);
        
        j_lfspd_t = -reshape(repmat(timeFacCurr',apertureInfo.beam(1).numOfActiveLeafPairs,2).*abs([leftLeafDiff rightLeafDiff])./ ...
            repmat((timeDoseBorderAngles.^2)',apertureInfo.beam(1).numOfActiveLeafPairs,2),2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeDoseBorderAngles),1);
        
        j_lfspd_cur(leftLeafDiff == 0) = abs(j_lfspd_cur(leftLeafDiff == 0));
        j_lfspd_nxt(leftLeafDiff == 0) = abs(j_lfspd_nxt(leftLeafDiff == 0));
        j_lfspd_cur(rightLeafDiff == 0) = abs(j_lfspd_cur(rightLeafDiff == 0));
        j_lfspd_nxt(rightLeafDiff == 0) = abs(j_lfspd_nxt(rightLeafDiff == 0));
        
        s = [j_lfspd_cur; j_lfspd_nxt; j_lfspd_t];
        
        jacob_lfspd = sparse(i,j,s,2*apertureInfo.beam(1).numOfActiveLeafPairs*apertureInfo.totalNumOfShapes,numel(apertureInfoVec),6*apertureInfo.beam(1).numOfActiveLeafPairs*apertureInfo.totalNumOfShapes);
        
        
        % jacobian of the doserate constraint
        % values of doserate (MU/sec) between optimized gantry angles
        weights = apertureInfoVec(1:apertureInfo.totalNumOfShapes)./apertureInfo.jacobiScale(1:apertureInfo.totalNumOfShapes);
        
        i = repmat(1:apertureInfo.totalNumOfShapes,1,2);
        j = [1:apertureInfo.totalNumOfShapes (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2+1):(apertureInfo.totalNumOfShapes*2+apertureInfo.totalNumOfLeafPairs*2)];
        % first do jacob wrt weights, then wrt times
        
        s = [apertureInfo.weightToMU./(timeDoseBorderAngles.*apertureInfo.jacobiScale(1:apertureInfo.totalNumOfShapes)); -apertureInfo.weightToMU.*weights.*timeFacCurr./(timeDoseBorderAngles.^2)];
        
        jacob_dosrt = sparse(i,j,s,apertureInfo.totalNumOfShapes,numel(apertureInfoVec),2*(apertureInfo.totalNumOfShapes));
    else
        
        % get index values for the jacobian
        % variable index
        % value of constraints for leaves
        leftLeafPos  = apertureInfoVec((1:apertureInfo.totalNumOfLeafPairs)+apertureInfo.totalNumOfShapes);
        rightLeafPos = apertureInfoVec(1+apertureInfo.totalNumOfLeafPairs+apertureInfo.totalNumOfShapes:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2);
        
        % values of time differences of optimized gantry angles
        DAOInd = [apertureInfo.propVMAT.beam.DAOBeam];
        timeDAOBorderAngles = apertureInfoVec((apertureInfo.numPhases*(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)+1):end);
        
        i = sort(repmat(1:(apertureInfo.totalNumOfShapes-1),1,2));
        j = sort(repmat(1:apertureInfo.totalNumOfShapes,1,2));
        j(1) = [];
        j(end) = [];
        
        timeFac = [apertureInfo.propVMAT.beam(DAOInd).timeFac]';
        timeFac(timeFac == 0) = [];
        
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
        
        jacob_lfspd = sparse(i,j,s,2*apertureInfo.beam(1).numOfActiveLeafPairs*(apertureInfo.totalNumOfShapes-1),numel(apertureInfoVec),numel(s));
        
        
        % jacobian of the doserate constraint
        % values of doserate (MU/sec) between optimized gantry angles
        weights = apertureInfoVec(1:apertureInfo.totalNumOfShapes)./apertureInfo.jacobiScale(1:apertureInfo.totalNumOfShapes);
        timeFacCurr = [apertureInfo.propVMAT.beam(DAOInd).timeFacCurr]';
        timeOptDoseBorderAngles = timeDAOBorderAngles.*timeFacCurr;
        
        i = repmat(1:apertureInfo.totalNumOfShapes,1,2);
        j = [1:apertureInfo.totalNumOfShapes (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2+1):(apertureInfo.totalNumOfShapes*2+apertureInfo.totalNumOfLeafPairs*2)];
        % first do jacob wrt weights, then wrt times
        
        s = [apertureInfo.weightToMU./(timeOptDoseBorderAngles.*apertureInfo.jacobiScale(1:apertureInfo.totalNumOfShapes)); -apertureInfo.weightToMU.*weights.*timeFacCurr./(timeOptDoseBorderAngles.^2)];
        
        jacob_dosrt = sparse(i,j,s,apertureInfo.totalNumOfShapes,numel(apertureInfoVec),2*(apertureInfo.totalNumOfShapes));
    end
    
    
    % concatenate
    jacob = [jacob_dao; jacob_lfspd; jacob_dosrt; jacob_dos];
end


