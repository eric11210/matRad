function jacob = matRad_jacobFuncWrapper(w,dij,cst,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: jacobian function for inverse planning supporting max dose
% constraint, min dose constraint, min mean dose constraint, max mean dose constraint,
% min EUD constraint, max EUD constraint, max DVH constraint, min DVH constraint
%
% call
%   jacob = matRad_jacobFunc(w,dij,cst,options)
%
% input
%   w:    bixel weight vector
%   dij:  dose influence matrix
%   cst:  matRad cst struct
%   options: option struct defining the type of optimization
%
% output
%   jacob: jacobian of constraint function
%
% References
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get current dose / effect / RBExDose vector
d = matRad_backProjection(w,dij,options);

% initialize projection matrices and id containers
physicalDoseProjection  = sparse([]);
mAlphaDoseProjection    = sparse([]);
mSqrtBetaDoseProjection = sparse([]);
voxelID                 = [];
constraintID            = 0;
scenID                  = [];
scenID2                 = [];

for scen = 1:options.numOfScenarios
    % compute objective function for every VOI.
    for i = 1:size(cst,1)
        
        % Only take OAR or target VOI.
        if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )
            
            % loop over the number of constraints for the current VOI
            for j = 1:numel(cst{i,6})
                
                % only perform computations for constraints
                if ~isempty(strfind(cst{i,6}(j).type,'constraint'))
                    
                    % compute reference
                    if (~isequal(cst{i,6}(j).type, 'max dose constraint')      && ~isequal(cst{i,6}(j).type, 'min dose constraint')          &&...
                            ~isequal(cst{i,6}(j).type, 'max mean dose constraint') && ~isequal(cst{i,6}(j).type, 'min mean dose constraint') && ...
                            ~isequal(cst{i,6}(j).type, 'min EUD constraint')       && ~isequal(cst{i,6}(j).type, 'max EUD constraint'))           && ...
                            isequal(options.bioOpt,'LEMIV_effect')
                        
                        d_ref = cst{i,5}.alphaX*cst{i,6}(j).dose + cst{i,5}.betaX*cst{i,6}(j).dose^2;
                    else
                        d_ref = cst{i,6}(j).dose;
                    end
                    
                    % if conventional opt: just add constraints of nominal dose
                    if strcmp(cst{i,6}(j).robustness,'none')
                        
                        d_i = d(cst{i,4}{1});
                        
                        jacobVec =  matRad_jacobFunc(d_i,cst{i,6}(j),d_ref);
                        
                        scenID  = [scenID;scen];
                        scenID2 = [scenID2;ones(numel(cst{i,4}{1}),1)];
                        
                        if isequal(options.bioOpt,'none') && ~isempty(jacobVec) || isequal(options.ID,'protons_const_RBExD')
                            
                            physicalDoseProjection = [physicalDoseProjection,sparse(cst{i,4}{1},1,jacobVec,dij.numOfVoxels,1)];
                            
                        elseif isequal(options.bioOpt,'LEMIV_effect') && ~isempty(jacobVec)
                            
                            mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{i,4}{1},1,jacobVec,dij.numOfVoxels,1)];
                            mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                sparse(cst{i,4}{1},1:numel(cst{i,4}{1}),2*jacobVec,dij.numOfVoxels,numel(cst{i,4}{1}))];
                            voxelID                 = [voxelID ;cst{i,4}{1}];
                            constraintID            = [constraintID, repmat(1 + constraintID(end),1,numel(cst{i,4}{1}))];
                            
                        elseif isequal(options.bioOpt,'LEMIV_RBExD') && ~isempty(jacobVec)
                            
                            scaledEffect = (dij.gamma(cst{i,4}{1}) + d_i);
                            
                            delta = jacobVec./(2*dij.bx(cst{i,4}{1}).*scaledEffect);
                            
                            mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{i,4}{1},1,delta,dij.numOfVoxels,1)];
                            mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                sparse(cst{i,4}{1},1:numel(cst{i,4}{1}),2*delta,dij.numOfVoxels,numel(cst{i,4}{1}))];
                            voxelID                 = [voxelID ;cst{i,4}{1}];
                            constraintID            = [constraintID, repmat(1 + constraintID(end),1,numel(cst{i,4}{1}))];
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

if isequal(options.bioOpt,'LEMIV_effect') || isequal(options.bioOpt,'LEMIV_RBExD')
    
    % initialize jacobian
    jacob = sparse([]);
    
    constraintID = constraintID(2:end);
end



numOfConstraints = numel(scenID)/options.numOfScenarios;
if isfield(dij,'optBixel')
    numOptBixel = nnz(dij.optBixel);
else
    numOptBixel = dij.totalNumOfBixels;
end

if options.FMO
    jacobSparseVec = zeros(1,numOfConstraints*numOptBixel*options.numOfScenarios);
else
    jacobSparseVec = cell(options.numOfScenarios,1);
    jacobSparseVec(:) = {zeros(1,numOfConstraints*numOptBixel)};
end

setOfConstraints = repmat(1:numOfConstraints,1,options.numOfScenarios);

% Calculate jacobian with dij projections
for i = 1:options.numOfScenarios
    % enter if statement also for protons using a constant RBE
    if isequal(options.bioOpt,'none') ||  isequal(options.ID,'protons_const_RBExD')
        
        if ~isempty(physicalDoseProjection)
            
            jacobLogical          = (scenID == i);
            currConstraints = setOfConstraints(jacobLogical);
            
            
            indInSparseVec = repmat(1:numOptBixel,1,numel(currConstraints))...
                +kron((currConstraints-1)*numOptBixel,ones(1,numOptBixel));
            
            offset = (i-1)*numOfConstraints*numOptBixel;
            
            if isfield(dij,'optBixel')
                
                if options.FMO
                    jacobSparseVec(offset+indInSparseVec) = transpose(physicalDoseProjection(:,jacobLogical)' * dij.scaleFactor * dij.physicalDose{i}(:,dij.optBixel));
                else
                    jacobSparseVec{i}(indInSparseVec) = transpose(physicalDoseProjection(:,jacobLogical)' * dij.scaleFactor * dij.physicalDose{i}(:,dij.optBixel));
                end
            else
                
                if options.FMO
                    jacobSparseVec(offset+indInSparseVec) = transpose(physicalDoseProjection(:,jacobLogical)' * dij.scaleFactor * dij.physicalDose{i});
                else
                    jacobSparseVec{i}(indInSparseVec) = transpose(physicalDoseProjection(:,jacobLogical)' * dij.scaleFactor * dij.physicalDose{i});
                end
            end
            
            if dij.memorySaverPhoton
                jacobVariables.currConstraints = currConstraints;
                jacobVariables.jacobLogical = jacobLogical;
                
                jacobSparseVec = jacobSparseVec+matRad_memorySaverDoseAndGrad(physicalDoseProjection,dij,'jacobian',i,jacobVariables);
            end
            
        end
        
    elseif isequal(options.bioOpt,'LEMIV_effect') || isequal(options.bioOpt,'LEMIV_RBExD')
        
        if ~isempty(mSqrtBetaDoseProjection) && ~isempty(mAlphaDoseProjection)
            
            jacobLogical            = (scenID == i);
            jacobLogical2           = (scenID2 == i);
            mSqrtBetaDoseProjection = mSqrtBetaDoseProjection(:,jacobLogical2)' * dij.mSqrtBetaDose{i} * w;
            mSqrtBetaDoseProjection = sparse(voxelID(jacobLogical2),constraintID(jacobLogical2),mSqrtBetaDoseProjection,...
                size(mAlphaDoseProjection(:,jacobLogical),1),size(mAlphaDoseProjection(:,jacobLogical),2));
            
            jacob(jacobLogical,:)   = mAlphaDoseProjection(:,jacobLogical)' * dij.mAlphaDose{i} +...
                mSqrtBetaDoseProjection' * dij.mSqrtBetaDose{i};
            
        end
    end
end

i_sparse = 1:numOfConstraints;
i_sparse = kron(i_sparse,ones(1,numOptBixel));

j_sparse = 1:dij.totalNumOfBixels;
if isfield(dij,'optBixel')
    j_sparse(~dij.optBixel) = [];
end
j_sparse = repmat(j_sparse,1,numOfConstraints);

if isequal(options.bioOpt,'none') ||  isequal(options.ID,'protons_const_RBExD')
    if options.FMO
        
        i_sparse = repmat(i_sparse,1,options.numOfScenarios);
        j_sparse = repmat(j_sparse,1,options.numOfScenarios);
        
        jOffset = repelem(((1:options.numOfScenarios)-1)*dij.totalNumOfBixels,numOptBixel);
        jOffset = repelem(jOffset,1,numOfConstraints);
        j_sparse = j_sparse+jOffset;
        
        jacob = sparse(i_sparse,j_sparse,jacobSparseVec,numOfConstraints,options.numOfScenarios*dij.totalNumOfBixels);
    else
        
        jacob = cell(options.numOfScenarios,1);
        for i = 1:options.numOfScenarios
            jacob{i} = sparse(i_sparse,j_sparse,jacobSparseVec{i},numOfConstraints,dij.totalNumOfBixels);
        end
    end
end

end
