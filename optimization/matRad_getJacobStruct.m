function jacobStruct = matRad_getJacobStruct(dij,cst,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: jacobian structure function for inverse planning supporting max dose
% constraint, min dose constraint, min mean dose constraint, max mean dose constraint,
% min EUD constraint, max EUD constraint, max DVH constraint, min DVH constraint
%
% call
%   jacobStruct = matRad_getJacobStruct(dij,cst)
%
% input
%   dij: dose influence matrix
%   cst: matRad cst struct
%
% output
%   jacobStruct: jacobian of constraint function
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


numOfConstraints = 0;
% check consitency of constraints
for i = 1:size(cst,1)
    for j = 1:numel(cst{i,6})
        if isequal(cst{i,6}(j).type, 'max mean dose constraint') || ...
                isequal(cst{i,6}(j).type, 'min mean dose constraint')
            
            % -> no other max, min mean dose constraint
            for k = j+1:numel(cst{i,6})
                if isequal(cst{i,6}(k).type, 'max mean dose constraint') || ...
                        isequal(cst{i,6}(k).type, 'min mean dose constraint')
                    error('Simultatenous definition of min mean and max mean dose constraint\n');
                end
            end
            
        elseif isequal(cst{i,6}(j).type, 'max EUD constraint') || ...
                isequal(cst{i,6}(j).type, 'min EUD constraint')
            
            % -> no other max, min mean dose constraint
            for k = j+1:numel(cst{i,6})
                if isequal(cst{i,6}(k).type, 'max EUD constraint') || ...
                        isequal(cst{i,6}(k).type, 'min EUD constraint')
                    error('Simultatenous definition of min EUD and max EUD constraint\n');
                end
            end
            
        elseif isequal(cst{i,6}(j).type, 'max DVH constraint') ||...
                isequal(cst{i,6}(j).type, 'min DVH constraint')
            
            % -> no other DVH constraint
            for k = j+1:numel(cst{i,6})
                if (isequal(cst{i,6}(k).type, 'max DVH constraint') && isequal(cst{i,6}(k).dose,cst{i,6}(j).dose)) || ...
                        (isequal(cst{i,6}(k).type, 'max DVH constraint') && isequal(cst{i,6}(k).volume,cst{i,6}(j).volume)) || ...
                        (isequal(cst{i,6}(k).type, 'min DVH constraint') && isequal(cst{i,6}(k).dose,cst{i,6}(j).dose)) || ...
                        (isequal(cst{i,6}(k).type, 'min DVH constraint') && isequal(cst{i,6}(k).volume,cst{i,6}(j).volume))
                    
                    error('Simultatenous definition of DVH constraint\n');
                end
            end
        end
        
        if strfind(cst{i,6}(j).type,'constraint')
            numOfConstraints = numOfConstraints+1;
        end
    end
end

% Initializes constraints
if isfield(dij,'optBixel')
    numOptBixel = nnz(dij.optBixel);
else
    numOptBixel = dij.totalNumOfBixels;
end

if options.FMO
    jacobStructVec = zeros(1,numOfConstraints*numOptBixel*options.numOfScenarios);
else
    jacobStructVec = cell(options.numOfScenarios,1);
    jacobStructVec{:} = {zeros(1,numOfConstraints*numOptBixel)};
end

for scen = 1:options.numOfScenarios
    
    offset = 0;
    
    % compute objective function for every VOI.
    for i = 1:size(cst,1)
        
        % Only take OAR or target VOI.
        if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )
            
            % loop over the number of constraints for the current VOI
            for j = 1:numel(cst{i,6})
                
                % only perform computations for constraints
                if ~isempty(strfind(cst{i,6}(j).type,'constraint'))
                    
                    % if conventional opt: just add constraints of nominal dose
                    if strcmp(cst{i,6}(j).robustness,'none')
                        
                        if isequal(cst{i,6}(j).type, 'max dose constraint') || ...
                                isequal(cst{i,6}(j).type, 'min dose constraint') || ...
                                isequal(cst{i,6}(j).type, 'max mean dose constraint') || ...
                                isequal(cst{i,6}(j).type, 'min mean dose constraint') || ...
                                isequal(cst{i,6}(j).type, 'max EUD constraint') || ...
                                isequal(cst{i,6}(j).type, 'min EUD constraint') || ...
                                isequal(cst{i,6}(j).type, 'max DVH constraint') || ...
                                isequal(cst{i,6}(j).type, 'min DVH constraint')
                            
                            FMOoffset = offset+(scen-1)*numOfConstraints*numOptBixel;
                            
                            if isfield(dij,'optBixel')
                                if options.FMO
                                    jacobStructVec(FMOoffset+(1:nnz(dij.optBixel))) = mean(dij.physicalDose{1}(cst{i,4}{1},dij.optBixel));
                                else
                                    jacobStructVec{scen}(offset+(1:nnz(dij.optBixel))) = mean(dij.physicalDose{1}(cst{i,4}{1},dij.optBixel));
                                end
                            else
                                if options.FMO
                                    jacobStructVec(FMOoffset+(1:numOptBixel)) = mean(dij.physicalDose{1}(cst{i,4}{1},:));
                                else
                                    jacobStructVec{scen}(offset+(1:numOptBixel)) = mean(dij.physicalDose{1}(cst{i,4}{1},:));
                                end
                            end
                            
                            if dij.memorySaverPhoton
                                jacobVariables.voxInStructure = cst{i,4}{1};
                                jacobVariables.numOfConstraints = numOfConstraints;
                                jacobVariables.offset = offset;
                                
                                jacobStructVec = jacobStructVec+matRad_memorySaverDoseAndGrad([],dij,'jacobianStruct',jacobVariables);
                            end
                            
                            offset = offset+numOptBixel;
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

i = 1:numOfConstraints;
i = kron(i,ones(1,numOptBixel));

j = 1:dij.totalNumOfBixels;
if isfield(dij,'optBixel')
    j(~dij.optBixel) = [];
end
j = repmat(j,1,numOfConstraints);

i = repmat(i,1,options.numOfScenarios);
j = repmat(j,1,options.numOfScenarios);

jOffset = repelem(((1:options.numOfScenarios)-1)*dij.totalNumOfBixels,numOptBixel);
jOffset = repelem(jOffset,1,numOfConstraints);
j = j+jOffset;

jacobStructVec(jacobStructVec ~= 0) = 1;

jacobStruct = sparse(i,j,jacobStructVec,numOfConstraints,options.numOfScenarios*dij.totalNumOfBixels);
