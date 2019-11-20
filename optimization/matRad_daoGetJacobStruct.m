function jacobStruct = matRad_daoGetJacobStruct(apertureInfo,dij,cst,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: get jacobian structure for direct aperture optimization
%
% call
%   jacobStruct = matRad_daoGetJacobStruct(apertureInfo,dij,cst)
%
% input
%   apertureInfo: aperture info struct
%   dij:          dose influence matrix
%   cst:          matRad cst struct
%
% output
%   jacobStruct: jacobian of constraint function
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

% jacobian structure of the dao constraints
% row indices
i = repmat(1:apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases,1,2);
% column indices
j = [(apertureInfo.totalNumOfShapes*apertureInfo.numPhases+1):(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs)*apertureInfo.numPhases ...
    ((apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs)*apertureInfo.numPhases+1):(apertureInfo.totalNumOfShapes+2*apertureInfo.totalNumOfLeafPairs)*apertureInfo.numPhases];

% -1 for left leaves, 1 for right leaves
s = ones(1,2*apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases);

jacobStruct_dao = sparse(i,j,s, ...
    apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases, ...
    numel(apertureInfo.apertureVector), ...
    2*apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases);

jacobStruct_dos_bixel = matRad_getJacobStruct(dij,cst,options);
% --> gives me a matrix with number of rows = num of constraints and tells
% me in th columns if a beamlet has an influence on this constraint

% for apertures I need to check if the very beam orientation of the aperture has a bixel
% that potentially influences the constraint

% for leaves I need to check if that particular leaf row has bixels that
% potentially influence the objective which works via apertureInfo.beam(i).bixelIndMap

% all stuff can be done per beam direction and then I use repmat to build
% up the big matrix

numOfConstraints = size(jacobStruct_dos_bixel{1},1);

i_sparse = 1:numOfConstraints;
i_sparse = kron(i_sparse,ones(1,numel(apertureInfo.apertureVector)));

j_sparse = 1:numel(apertureInfo.apertureVector);
j_sparse = repmat(j_sparse,1,numOfConstraints);

jacobStructSparseVec = zeros(numOfConstraints*numel(apertureInfo.apertureVector),1);

offset = 1;
if apertureInfo.propVMAT.continuousAperture
    repFactor = 2;
else
    repFactor = 1;
end

for phase = 1:apertureInfo.numPhases
    
    if ~isempty(jacobStruct_dos_bixel{phase})
        
        if apertureInfo.runVMAT
            
            DAOBeams = find([apertureInfo.propVMAT.beam.DAOBeam]);
            
            for i = 1:numel(apertureInfo.beam)
                
                % get used bixels in beam
                ixWeight = ~isnan(apertureInfo.beam(i).bixelIndMap);
                
                if apertureInfo.propVMAT.beam(i).DAOBeam
                    % DAO beam, don't worry about adding since this is just
                    % struct, i.e. we are only interested if the element is
                    % non-zero
                    
                    % first weight
                    jacobStructSparseVec(offset == j_sparse) = jacobStructSparseVec(offset == j_sparse)+sum(jacobStruct_dos_bixel{phase}(:,apertureInfo.beam(i).bixelIndMap(ixWeight)),2);
                    
                    % now leaf positions
                    for k = 1:apertureInfo.beam(i).numOfActiveLeafPairs
                        
                        ixLeaf  = ~isnan(apertureInfo.beam(i).bixelIndMap(k,:));
                        indInBixVec = apertureInfo.beam(i).bixelIndMap(k,ixLeaf);
                        
                        indInOptVec = apertureInfo.beam(i).shape{phase}(1).vectorOffset+k-1;
                        indInOptVec = repmat(indInOptVec,1,repFactor)+repelem([0 apertureInfo.numPhases*apertureInfo.totalNumOfLeafPairs],1,repFactor);
                        
                        indInSparseVec = repmat(indInOptVec,1,numOfConstraints)...
                            +kron((0:numOfConstraints-1)*numel(apertureInfo.apertureVector),ones(1,2*repFactor));
                        
                        jacobStructSparseVec(indInSparseVec) = jacobStructSparseVec(indInSparseVec)+repelem(sum(jacobStruct_dos_bixel{phase}(:,indInBixVec),2),2*repFactor,1);
                    end
                    
                    offset = offset+1;
                else
                    % not DAO beam, these may contain bixels which affect the
                    % constraints which are influenced by DAO leaf pairs that
                    % do not affect the constraints (unlikely to happen, but it
                    % might)
                    
                    %first weight
                    
                    %give fraction of gradient to previous optimized beam
                    lastDAOInd = find(DAOBeams == apertureInfo.propVMAT.beam(i).lastDAOIndex,1)+(phase-1)*apertureInfo.totalNumOfShapes;
                    jacobStructSparseVec(lastDAOInd == j_sparse) = jacobStructSparseVec(lastDAOInd == j_sparse)+sum(jacobStruct_dos_bixel{phase}(:,apertureInfo.beam(i).bixelIndMap(ixWeight)),2);
                    %give the other fraction to next optimized beam
                    nextDAOInd = find(DAOBeams == apertureInfo.propVMAT.beam(i).nextDAOIndex,1)+(phase-1)*apertureInfo.totalNumOfShapes;
                    jacobStructSparseVec(nextDAOInd == j_sparse) = jacobStructSparseVec(nextDAOInd == j_sparse)+sum(jacobStruct_dos_bixel{phase}(:,apertureInfo.beam(i).bixelIndMap(ixWeight)),2);
                    
                    %now leaf pos
                    
                    for k = 1:apertureInfo.beam(i).numOfActiveLeafPairs
                        
                        ixLeaf  = ~isnan(apertureInfo.beam(i).bixelIndMap(k,:));
                        indInBixVec = apertureInfo.beam(i).bixelIndMap(k,ixLeaf);
                        
                        %give fraction of gradient to previous optimized beam
                        indInOptVec = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape{phase}(1).vectorOffset+k-1;
                        indInOptVec = repmat(indInOptVec,1,repFactor)+repelem([0 apertureInfo.numPhases*apertureInfo.totalNumOfLeafPairs],1,repFactor);
                        indInSparseVec = repmat(indInOptVec,1,numOfConstraints)...
                            +kron((0:numOfConstraints-1)*numel(apertureInfo.apertureVector),ones(1,2*repFactor));
                        
                        jacobStructSparseVec(indInSparseVec) = jacobStructSparseVec(indInSparseVec)+repelem(sum(jacobStruct_dos_bixel{phase}(:,indInBixVec),2),2*repFactor,1);
                        
                        
                        %give the other fraction to next optimized beam
                        indInOptVec = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape{phase}(1).vectorOffset+k-1;
                        indInOptVec = repmat(indInOptVec,1,repFactor)+repelem([0 apertureInfo.numPhases*apertureInfo.totalNumOfLeafPairs],1,repFactor);
                        indInSparseVec = repmat(indInOptVec,1,numOfConstraints)...
                            +kron((0:numOfConstraints-1)*numel(apertureInfo.apertureVector),ones(1,2*repFactor));
                        
                        jacobStructSparseVec(indInSparseVec) = jacobStructSparseVec(indInSparseVec)+repelem(sum(jacobStruct_dos_bixel{phase}(:,indInBixVec),2),2*repFactor,1);
                    end
                    
                    %now time
                    
                    %give fraction of gradient to previous optimized beam
                    lastDAOIndTime = find(DAOBeams == apertureInfo.propVMAT.beam(i).lastDAOIndex,1)+(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases;
                    jacobStructSparseVec(lastDAOIndTime == j_sparse) = jacobStructSparseVec(lastDAOIndTime == j_sparse)+sum(jacobStruct_dos_bixel{phase}(:,apertureInfo.beam(i).bixelIndMap(ixWeight)),2);
                    
                    %give the other fraction to next optimized beam
                    nextDAOIndTime = find(DAOBeams == apertureInfo.propVMAT.beam(i).nextDAOIndex,1)+(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases;
                    jacobStructSparseVec(nextDAOIndTime == j_sparse) = jacobStructSparseVec(nextDAOIndTime == j_sparse)+sum(jacobStruct_dos_bixel{phase}(:,apertureInfo.beam(i).bixelIndMap(ixWeight)),2);
                end
            end
        else
            
            for i = 1:numel(apertureInfo.beam)
                
                % get used bixels in beam
                ixWeight = ~isnan(apertureInfo.beam(i).bixelIndMap);
                
                for j = 1:apertureInfo.beam(i).numOfShapes
                    % first weight
                    jacobStructSparseVec(offset+j == j_sparse) = jacobStructSparseVec(offset+j == j_sparse)+sum(jacobStruct_dos_bixel(:,apertureInfo.beam(i).bixelIndMap(ixWeight)),2);
                    
                    % now leaf positions
                    for k = 1:apertureInfo.beam(i).numOfActiveLeafPairs
                        
                        ixLeaf  = ~isnan(apertureInfo.beam(i).bixelIndMap(k,:));
                        indInBixVec = apertureInfo.beam(i).bixelIndMap(k,ixLeaf);
                        
                        indInOptVec = apertureInfo.beam(i).shape(1).vectorOffset+k-1+[0 apertureInfo.totalNumOfLeafPairs];
                        indInSparseVec = repmat(indInOptVec,1,numOfConstraints)...
                            +kron((0:numOfConstraints-1)*numel(apertureInfo.apertureVector),ones(1,2));
                        
                        jacobStructSparseVec(indInSparseVec) = jacobStructSparseVec(indInSparseVec)+repelem(sum(jacobStruct_dos_bixel(:,indInBixVec),2),2,1);
                    end
                    
                    offset = offset+1;
                end
            end
        end
        
        jacobStructSparseVec(jacobStructSparseVec ~= 0) = 1;
        
        jacobStruct_dos = sparse(i_sparse,j_sparse,jacobStructSparseVec,numOfConstraints,numel(apertureInfo.apertureVector));
    else
        jacobStruct_dos = sparse(0,0);
    end
    
end


if ~apertureInfo.runVMAT
    % concatenate
    jacobStruct = [jacobStruct_dao; jacobStruct_dos];
else
    
    if apertureInfo.propVMAT.continuousAperture
        % set up
        n = apertureInfo.beam(1).numOfActiveLeafPairs;
        indInSparseVec  = (1:n);
        indInConVec     = (1:n);
        shapeInd        = 1;
        
        % sparse matrix
        numElem     = n.*apertureInfo.propVMAT.numLeafSpeedConstraintDAO*6;
        i_sparse    = zeros(numElem,1);
        j_sparse    = zeros(numElem,1);
        s_sparse    = ones(numElem,1);
        
        for i = 1:numel(apertureInfo.beam)
            % loop over beams
            
            if ~isempty(apertureInfo.propVMAT.beam(i).leafConstMask)
                % only true for DAO beams
                for phase_I = 1:apertureInfo.numPhases
                    % loop over initial phases
                    
                    transitions = apertureInfo.propVMAT.beam(i).leafConstMask(phase_I,:);
                    transitions(transitions == 0) = [];
                    
                    for phase_F = transitions
                        % loop over possible final phases
                        
                        % get vector indices
                        % it's a DAO beam, so use own vector offset
                        vectorIx_LI = apertureInfo.beam(i).shape{phase_I}(1).vectorOffset(1) + ((1:n)-1);
                        vectorIx_LF = apertureInfo.beam(i).shape{phase_F}(1).vectorOffset(2) + ((1:n)-1);
                        vectorIx_RI = vectorIx_LI+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
                        vectorIx_RF = vectorIx_LF+apertureInfo.totalNumOfLeafPairs*apertureInfo.numPhases;
                        
                        % calc jacobs
                        
                        % wrt initial leaf positions (left, then right)
                        i_sparse(indInSparseVec)    = indInConVec;
                        j_sparse(indInSparseVec)    = vectorIx_LI;
                        indInSparseVec              = indInSparseVec+n;
                        
                        i_sparse(indInSparseVec)    = indInConVec+apertureInfo.propVMAT.numLeafSpeedConstraint*apertureInfo.beam(1).numOfActiveLeafPairs;
                        j_sparse(indInSparseVec)    = vectorIx_RI;
                        indInSparseVec              = indInSparseVec+n;
                        
                        % wrt final leaf positions (left, then right)
                        i_sparse(indInSparseVec)    = indInConVec;
                        j_sparse(indInSparseVec)    = vectorIx_LF;
                        indInSparseVec              = indInSparseVec+n;
                        
                        i_sparse(indInSparseVec)    = indInConVec+apertureInfo.propVMAT.numLeafSpeedConstraint*apertureInfo.beam(1).numOfActiveLeafPairs;
                        j_sparse(indInSparseVec)    = vectorIx_RF;
                        indInSparseVec              = indInSparseVec+n;
                        
                        % wrt time (left, then right)
                        % it's a DAO beam, so speed only depends on its own
                        % time
                        i_sparse(indInSparseVec)    = indInConVec;
                        j_sparse(indInSparseVec)    = apertureInfo.propVMAT.beam(i).timeInd;
                        indInSparseVec              = indInSparseVec+n;
                        
                        i_sparse(indInSparseVec)    = indInConVec+apertureInfo.propVMAT.numLeafSpeedConstraint*apertureInfo.beam(1).numOfActiveLeafPairs;
                        j_sparse(indInSparseVec)    = apertureInfo.propVMAT.beam(i).timeInd;
                        indInSparseVec              = indInSparseVec+n;
                        
                        % update offset
                        indInConVec = indInConVec+n;
                    end
                end
                
                % increment shapeInd only for beams which have transtion
                % defined
                shapeInd = shapeInd+1;
            end
        end
        
        jacobStruct_lfspd = sparse(i_sparse,j_sparse,s_sparse,2*apertureInfo.beam(1).numOfActiveLeafPairs*apertureInfo.propVMAT.numLeafSpeedConstraint,numel(apertureInfo.apertureVector));
        
        %{
        % get index values for the jacobian
        % variable index
        initialLeftLeafInd = [apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).initialLeftLeafInd];
        initialRightLeafInd = [apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).initialRightLeafInd];
        finalLeftLeafInd = [apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).finalLeftLeafInd];
        finalRightLeafInd = [apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).finalRightLeafInd];
        leftTimeInd = [apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).timeInd];
        rightTimeInd = [apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).timeInd];
        % constraint index
        leafConstraintInd = 1:2*apertureInfo.beam(1).numOfActiveLeafPairs*apertureInfo.propVMAT.numLeafSpeedConstraint;
        
        timeFacInd = [apertureInfo.propVMAT.beam.timeFacInd];
        timeFacInd([apertureInfo.propVMAT.beam.timeFac] == 0) = [];
        
        timeConstraintInd = repmat(1:apertureInfo.beam(1).numOfActiveLeafPairs,1,apertureInfo.propVMAT.numLeafSpeedTimeEffect)+repelem((timeFacInd-1)*apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.beam(1).numOfActiveLeafPairs);
        timeConstraintInd = repmat(timeConstraintInd,1,2)+repelem([0 apertureInfo.beam(1).numOfActiveLeafPairs*apertureInfo.propVMAT.numLeafSpeedConstraint],1,apertureInfo.beam(1).numOfActiveLeafPairs*apertureInfo.propVMAT.numLeafSpeedTimeEffect);
        
        % jacobian of the leafspeed constraint
        i = [repmat(leafConstraintInd,1,2) timeConstraintInd];
        j = [initialLeftLeafInd initialRightLeafInd finalLeftLeafInd finalRightLeafInd leftTimeInd rightTimeInd];
        % first do jacob wrt initial leaf position (left, right), then final leaf
        % position (left, right), then time (left, right)
        
        s = ones(1,(4*apertureInfo.propVMAT.numLeafSpeedConstraint+2*apertureInfo.propVMAT.numLeafSpeedTimeEffect)*apertureInfo.beam(1).numOfActiveLeafPairs);
        
        jacobStruct_lfspd = zeros(2*apertureInfo.propVMAT.numLeafSpeedConstraint*apertureInfo.beam(1).numOfActiveLeafPairs,numel(apertureInfo.apertureVector));
        %jacobStruct_lfspd = sparse(i,j,s,2*apertureInfo.beam(1).numOfActiveLeafPairs*apertureInfo.propVMAT.numLeafSpeedConstraint,numel(apertureInfo.apertureVector));
        %}
    else
        
        i = sort(repmat(1:(apertureInfo.totalNumOfShapes-1),1,2));
        j = sort(repmat(1:apertureInfo.totalNumOfShapes,1,2));
        j(1) = [];
        j(end) = [];
        
        % get index values for the jacobian
        % variable index
        timeInd = (1+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2):(apertureInfo.totalNumOfShapes*2+apertureInfo.totalNumOfLeafPairs*2-1);
        currentLeftLeafInd = (apertureInfo.totalNumOfShapes+1):(apertureInfo.totalNumOfShapes+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeInd));
        currentRightLeafInd = (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+1):(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeInd));
        nextLeftLeafInd = (apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+1):(apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeInd));
        nextRightLeafInd = (apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+1):(apertureInfo.beam(1).numOfActiveLeafPairs+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeInd));
        leftTimeInd = kron(j,ones(1,apertureInfo.beam(1).numOfActiveLeafPairs))+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2;
        rightTimeInd = kron(j,ones(1,apertureInfo.beam(1).numOfActiveLeafPairs))+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2;
        %leftTimeInd = repelem(j,apertureInfo.beam(1).numOfActiveLeafPairs)+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2;
        %rightTimeInd = repelem(j,apertureInfo.beam(1).numOfActiveLeafPairs)+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2;
        % constraint index
        leafConstraintInd = 1:2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeInd);
        
        % jacobian of the leafspeed constraint
        i = repmat((i'-1)*apertureInfo.beam(1).numOfActiveLeafPairs,1,apertureInfo.beam(1).numOfActiveLeafPairs)+repmat(1:apertureInfo.beam(1).numOfActiveLeafPairs,2*numel(timeInd),1);
        i = reshape([i' i'+apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeInd)],1,[]);
        
        i = [repmat(leafConstraintInd,1,2) i];
        j = [currentLeftLeafInd currentRightLeafInd nextLeftLeafInd nextRightLeafInd leftTimeInd rightTimeInd];
        % first do jacob wrt current leaf position (left, right), then next leaf
        % position (left, right), then time (left, right)
        
        s = ones(1,numel(j));
        
        jacobStruct_lfspd = sparse(i,j,s,2*apertureInfo.beam(1).numOfActiveLeafPairs*(apertureInfo.totalNumOfShapes-1),numel(apertureInfo.apertureVector),numel(s));
    end
    
    % jacobian of the doserate constraint
    i = repmat(1:(apertureInfo.totalNumOfShapes*apertureInfo.numPhases),1,2);
    j = [1:(apertureInfo.totalNumOfShapes*apertureInfo.numPhases) ...
        repmat(((apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases+1):((apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases+apertureInfo.totalNumOfShapes),1,apertureInfo.numPhases)];
    % first do jacob wrt weights, then wrt times
    
    s = ones(1,2*apertureInfo.totalNumOfShapes*apertureInfo.numPhases);
    
    jacobStruct_dosrt = sparse(i,j,s,apertureInfo.totalNumOfShapes*apertureInfo.numPhases,numel(apertureInfo.apertureVector),2*apertureInfo.totalNumOfShapes*apertureInfo.numPhases);
    
    
    % concatenate
    jacobStruct = [jacobStruct_dao; jacobStruct_lfspd; jacobStruct_dosrt; jacobStruct_dos];
end


