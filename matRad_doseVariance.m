function [dVarSum,dVarSumGrad] = matRad_doseVariance(apertureInfo,dij)
%#codegen

%% setup

% allocate raw doses
dRaw_I      = cell(apertureInfo.numPhases*numel(apertureInfo.beam),1);
dRaw_F      = cell(apertureInfo.numPhases*numel(apertureInfo.beam),1);

% allocate gradient of raw doses
dRawGradCell_I  = cell(apertureInfo.numPhases*numel(apertureInfo.beam),1);
dRawGradCell_F  = cell(apertureInfo.numPhases*numel(apertureInfo.beam),1);

% allocate mean of dose
d = zeros(numel(dij.targetVox),1);

% allocate mean of squared dose
d2Sum = 0;

% allocate gradient of mean of dose
dGrad = zeros(numel(dij.targetVox),size(apertureInfo.apertureVector,1));

% allocate gradient of mean of squared dose
d2SumGrad = zeros(1,size(apertureInfo.apertureVector,1));

% allocate time gradients of mean of squared dose
d2SumTimeGrad = zeros(1,apertureInfo.totalNumOfShapes);

% determine time derivative variable indices
tIxBig_Vec  = (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases+(1:apertureInfo.totalNumOfShapes);
tIx_Vec     = 1:apertureInfo.totalNumOfShapes;

%% precalculate doses and gradients
% loop over all beams

% track nnz
nnzGrad_I = 0;
nnzGrad_F = 0;

for i = 1:numel(apertureInfo.beam)
    
    % pre-calculate raw doses
    % first find relevant bixels
    currBixelIx = apertureInfo.beam(i).bixelIndMap(~isnan(apertureInfo.beam(i).bixelIndMap));
    
    % find relevant variables for beam
    currVarIx = apertureInfo.beam(i).local2GlobalVar;
    
    % allocate gradient of dose for beam
    dGrad_i = zeros(numel(dij.targetVox),apertureInfo.beam(i).numUniqueVar);
    
    % loop over phases
    for phase = 1:apertureInfo.numPhases
        
        % extract dijBeamPhase
        dijBeamPhase    = dij.scaleFactor .* dij.physicalDose{phase}(dij.targetVox,currBixelIx);
        dijBeamPhaseInv = dijBeamPhase';
        
        % now extract bixel weights
        w_I = apertureInfo.arcI.bixelWeights{phase}(currBixelIx)';
        w_F = apertureInfo.arcF.bixelWeights{phase}(currBixelIx)';
        % calculate dose
        dRaw_I{(i-1).*apertureInfo.numPhases+phase} = w_I * dijBeamPhaseInv;
        dRaw_F{(i-1).*apertureInfo.numPhases+phase} = w_F * dijBeamPhaseInv;
        
        % calculate gradient of dose
        % d = d + dij.scaleFactor .*
        % dij.physicalDose{phase_I1}(:,currBixelIx1) * apertureInfo.bixelWeights{phase_I1}(currBixelIx1);
        j_I = apertureInfo.arcI.bixelJApVec{(i-1).*apertureInfo.numPhases+phase}';
        j_F = apertureInfo.arcF.bixelJApVec{(i-1).*apertureInfo.numPhases+phase}';
        % calculate gradient
        dRawGradCell_I{(i-1).*apertureInfo.numPhases+phase} = dijBeamPhase * j_I;
        dRawGradCell_F{(i-1).*apertureInfo.numPhases+phase} = dijBeamPhase * j_F;
        
        % USE REGULAR ARRAY FOR dRawGrad INSTEAD OF CELL?
        % sum([apertureInfo.beam.numUniqueVar])*apertureInfo.numPhases
        
        % track number of nonzero elements
        nnzGrad_I = nnzGrad_I+nnz(dRawGradCell_I{(i-1).*apertureInfo.numPhases+phase});
        nnzGrad_F = nnzGrad_F+nnz(dRawGradCell_F{(i-1).*apertureInfo.numPhases+phase});
        
        % find relevant variables for gradient
        currVarIx = apertureInfo.beam(i).local2GlobalVar;
        
        % accumulate sum of dose and dose gradient
        d = d+dRaw_I{(i-1).*apertureInfo.numPhases+phase}'+dRaw_F{(i-1).*apertureInfo.numPhases+phase}';
        dGrad_i = dGrad_i+dRawGradCell_I{(i-1).*apertureInfo.numPhases+phase}+dRawGradCell_F{(i-1).*apertureInfo.numPhases+phase};
        
    end
    
    % dump gradient for i
    dGrad(:,currVarIx) = dGrad(:,currVarIx)+dGrad_i;
    
end

% preallocate full values and indices
V_I = zeros(nnzGrad_I,1);
I_I = zeros(nnzGrad_I,1);
J_I = zeros(nnzGrad_I,1);
V_F = zeros(nnzGrad_F,1);
I_F = zeros(nnzGrad_F,1);
J_F = zeros(nnzGrad_F,1);
%dRawGrad_I = spalloc(apertureInfo.beam(end).gradOffset(end)+apertureInfo.beam(end).numUniqueVar-1,numel(dij.targetVox),nnzGrad_I);
%dRawGrad_F = spalloc(apertureInfo.beam(end).gradOffset(end)+apertureInfo.beam(end).numUniqueVar-1,numel(dij.targetVox),nnzGrad_F);

% initialize offset
indOffset_I = 0;
indOffset_F = 0;

% loop over all beams again to put individual values and indices into full
for i = 1:numel(apertureInfo.beam)
    
    % loop over phases
    for phase = 1:apertureInfo.numPhases
        
        % extract gradients
        dRawGradTemp_I = dRawGradCell_I{(i-1).*apertureInfo.numPhases+phase};
        dRawGradTemp_F = dRawGradCell_F{(i-1).*apertureInfo.numPhases+phase};
        
        % get nonzero values and indices
        v_I = nonzeros(dRawGradTemp_I);
        v_F = nonzeros(dRawGradTemp_F);
        [i_I,j_I] = find(dRawGradTemp_I);
        [i_F,j_F] = find(dRawGradTemp_F);
        
        % get indices for full
        ind_I = indOffset_I+(1:numel(v_I));
        ind_F = indOffset_F+(1:numel(v_F));
        
        % put values and indives into full
        V_I(ind_I) = v_I;
        I_I(ind_I) = i_I;
        J_I(ind_I) = j_I;
        V_F(ind_F) = v_F;
        I_F(ind_F) = i_F;
        J_F(ind_F) = j_F;
        
        % wipe cell
        dRawGradCell_I{(i-1).*apertureInfo.numPhases+phase} = [];
        dRawGradCell_F{(i-1).*apertureInfo.numPhases+phase} = [];
        
        % increase offset
        indOffset_I = indOffset_I+numel(v_I);
        indOffset_F = indOffset_F+numel(v_F);
        
    end
end

% TRY SPARSE MATRIX AGAIN, THIS TIME WITH THE VOXELS AS THE ROW INDEX AND
% THE VARIABLES AS THE COLUMNS

% construct sparse matrix of gradients
dRawGrad_I = sparse(I_I,J_I,V_I,numel(dij.targetVox),apertureInfo.beam(end).gradOffset+apertureInfo.numPhases.*apertureInfo.beam(end).numUniqueVar,nnzGrad_I);
dRawGrad_F = sparse(I_F,J_F,V_F,numel(dij.targetVox),apertureInfo.beam(end).gradOffset+apertureInfo.numPhases.*apertureInfo.beam(end).numUniqueVar,nnzGrad_F);

% clear unnecessary variables
clear I_I J_I V_I I_F J_F V_F

% first loop over all beams
for i1 = 1:numel(apertureInfo.beam)
    
    % find relevant variables for beam 1
    currVarIx1 = apertureInfo.beam(i1).local2GlobalVar;
    
    % find kept variables for beam 1
    d2KeepVar1 = apertureInfo.beam(i1).d2KeepVar;
    %d2KeepVar1(:) = true;
    %d2KeepVar1(6:end) = false;
    
    % find variable indices for beam 2
    numUniqueVar_i1 = apertureInfo.beam(i1).numUniqueVar;
    varInd_i1 = apertureInfo.beam(i1).gradOffset+(1:apertureInfo.beam(i1).numUniqueVar);
    varInd_i1(~d2KeepVar1) = [];
    
    % allocate gradient of mean of squared dose for beam 1
    d2SumGrad_i1 = zeros(1,nnz(d2KeepVar1));
    %d2SumGrad_i1 = zeros(1,apertureInfo.beam(i1).numUniqueVar);
    
    % second loop over all beams
    for i2 = 1:numel(apertureInfo.beam)
        
        %% determine probabilities
        % 1) probability to transfer from the end of the earlier to the
        % beginning of the later
        % 2) probability to arrive at the later
        
        % determine which i is first and last
        iFirst = min([i1 i2]);
        iLast = max([i1 i2]);
        
        % determine pNormGradFactors
        pNormGradFactor1 = sum(apertureInfo.propVMAT.jacobT(:,(iFirst+1):(iLast-1)),2);
        pNormGradFactor2 = sum(apertureInfo.propVMAT.jacobT(:,1:(iLast-1)),2);
        
        % calculate the transition and arrival times
        transT12    = sum([apertureInfo.beam((iFirst+1):(iLast-1)).time]);
        T2          = sum([apertureInfo.beam(1:(iLast-1)).time]);
        
        % calculate the probabilities
        [Pij_transT12,Pij_transT12_dot,Pi_T2,Pi_T2_dot] = matRad_transAndTProb(apertureInfo.propVMAT.qij,apertureInfo.propVMAT.initProb,transT12,T2);
        
        % find relevant variables for beam 2
        currVarIx2 = apertureInfo.beam(i2).local2GlobalVar;
        
        % find kept variables for beam 2
        d2KeepVar2 = apertureInfo.beam(i2).d2KeepVar;
        %d2KeepVar2(:) = true;
        %d2KeepVar2(6:end) = false;
        
        % find variable indices for beam 2
        numUniqueVar_i2 = apertureInfo.beam(i2).numUniqueVar;
        varInd_i2 = apertureInfo.beam(i2).gradOffset+(1:apertureInfo.beam(i2).numUniqueVar);
        varInd_i2(~d2KeepVar2) = [];
        
        % allocate gradient of mean of squared dose for each beam
        %d2Grad_i1 = zeros(numel(dij.targetVox),apertureInfo.beam(i1).numUniqueVar);
        %d2Grad_i2 = zeros(numel(dij.targetVox),apertureInfo.beam(i2).numUniqueVar);
        
        % allocate gradient of mean of squared dose for beam 2
        d2SumGrad_i2 = zeros(1,nnz(d2KeepVar2));
        %d2SumGrad_i2 = zeros(1,apertureInfo.beam(i2).numUniqueVar);
        
        % first loop over initial phases
        for phase_I1 = 1:apertureInfo.numPhases
            
            % get pre-calculated raw doses
            dRaw_I1 = dRaw_I{(i1-1).*apertureInfo.numPhases+phase_I1};
            
            % get pre-calculated raw gradients
            % only consider gradients wrt weights and times for d2
            varInd              = (phase_I1-1).*numUniqueVar_i1+varInd_i1;
            dRawGrad_I1         = dRawGrad_I(:,varInd);
            %dRawGrad_I1 = dRawGradCell_I{(i1-1).*apertureInfo.numPhases+phase_I1}(:,d2KeepVar1);
            
            % first loop over final phases
            for phase_F1 = 1:apertureInfo.numPhases
                
                % get pre-calculated raw doses
                dRaw_F1 = dRaw_F{(i1-1).*apertureInfo.numPhases+phase_F1};
                
                % get pre-calculated raw gradients
                % only consider gradients wrt weights and times for d2
                varInd              = (phase_F1-1).*numUniqueVar_i1+varInd_i1;
                dRawGrad_F1         = dRawGrad_F(:,varInd);
                %dRawGrad_F1 = dRawGradCell_F{(i1-1).*apertureInfo.numPhases+phase_F1}(:,d2KeepVar1);
                
                % sum I1 and F1 doses and gradients
                dRaw_1      = dRaw_I1+dRaw_F1;
                dRawGrad_1  = dRawGrad_I1+dRawGrad_F1;
                
                % second loop over initial phases
                for phase_I2 = 1:apertureInfo.numPhases
                    
                    % determine probability normalization
                    pNorm = Pij_transT12(phase_F1,phase_I2)./Pi_T2(phase_I2);
                    
                    % calculate gradient of normalizing probability wrt
                    % times, update gradient vector
                    pNormGrad = (pNormGradFactor1.*Pij_transT12_dot(phase_F1,phase_I2)./Pi_T2(phase_I2) - ...
                        pNormGradFactor2.*Pi_T2_dot(phase_I2).*Pij_transT12(phase_F1,phase_I2)./(Pi_T2(phase_I2).^2))';
                    
                    % determine if we're doing d2 for this combo
                    dod2 = abs(pNorm) > eps;
                    
                    % determine if we're doing time gradients of d2 for
                    % this combo
                    deleteInd                   = abs(pNormGrad) < eps;
                    pNormGrad(deleteInd)        = [];
                    temptIx_Vec                 = tIx_Vec;
                    temptIx_Vec(deleteInd)      = [];
                    dod2TimeGrad                = ~all(deleteInd);
                    
                    % these are only necessary if we are doing d2 or
                    % d2TimeGrad stuff
                    if dod2 || dod2TimeGrad
                        % get pre-calculated raw doses
                        dRaw_I2 = dRaw_I{(i2-1).*apertureInfo.numPhases+phase_I2};
                        
                        if dod2
                            % get pre-calculated raw gradients
                            % only consider gradients wrt weights and times for d2
                            varInd              = (phase_I2-1).*numUniqueVar_i2+varInd_i2;
                            dRawGrad_I2         = dRawGrad_I(:,varInd);
                            %dRawGrad_I2 = dRawGradCell_I{(i2-1).*apertureInfo.numPhases+phase_I2}(:,d2KeepVar2);
                        end
                    end
                    
                    % second loop over final phases
                    for phase_F2 = 1:apertureInfo.numPhases
                        
                        %% determine mean of squared dose for this combination
                        
                        % check if we are doing either d2 or dod2TimeGrad stuff
                        if dod2 || dod2TimeGrad
                            
                            % get pre-calculated raw doses
                            dRaw_F2 = dRaw_F{(i2-1).*apertureInfo.numPhases+phase_F2};
                            
                            % sum I2 and F2 doses
                            dRaw_2 = dRaw_I2+dRaw_F2;
                            
                            % the "weight" tensors are given by of all
                            % combinations of weight products multiplied by the
                            % extra probability to observe both trajectories
                            % SKIP tensor determination step and go right to
                            % calculating mean squared doses
                            %d2SumRaw_II = dRaw_I1 * dRaw_I2';
                            %d2SumRaw_IF = dRaw_I1 * dRaw_F2';
                            %d2SumRaw_FI = dRaw_F1 * dRaw_I2';
                            %d2SumRaw_FF = dRaw_F1 * dRaw_F2';
                            
                            % combine into a single multiplication
                            d2SumRaw = dRaw_1 * dRaw_2';
                            
                            % only do the following if we need d2 stuff,
                            % i.e. if pNorm is > eps
                            if dod2
                                
                                % the "weight" tensors are given by of all
                                % combinations of weight products multiplied by the
                                % extra probability to observe both trajectories
                                % SKIP tensor determination step and go right to
                                % calculating mean squared doses
                                %d2Sum_II = pNorm .* d2SumRaw_II;
                                %d2Sum_IF = pNorm .* d2SumRaw_IF;
                                %d2Sum_FI = pNorm .* d2SumRaw_FI;
                                %d2Sum_FF = pNorm .* d2SumRaw_FF;
                                
                                % sum terms
                                %d2Sum  = d2Sum+d2Sum_II+d2Sum_IF+d2Sum_FI+d2Sum_FF;
                                d2Sum   = d2Sum + pNorm .* d2SumRaw;
                                
                                % get pre-calculated raw gradients
                                % only consider gradients wrt weights and times for d2
                                varInd              = (phase_F2-1).*numUniqueVar_i2+varInd_i2;
                                dRawGrad_F2         = dRawGrad_F(:,varInd);
                                %dRawGrad_F2 = dRawGradCell_F{(i2-1).*apertureInfo.numPhases+phase_F2}(:,d2KeepVar2);
                                
                                % sum I2 and F2 gradients
                                dRawGrad_2  = dRawGrad_I2+dRawGrad_F2;
                                
                                % calculate gradients of mean squared doses, put
                                % them directly in d2Grad
                                % first only look at gradients of the raw
                                % doses, not pNorm
                                
                                % do II
                                %d2SumGrad_i1   = d2SumGrad_i1 + pNorm .* dRaw_I2 * dRawGrad_I1;
                                %d2SumGrad_i2   = d2SumGrad_i2 + pNorm .* dRaw_I1 * dRawGrad_I2;
                                % do IF
                                %d2SumGrad_i1   = d2SumGrad_i1 + pNorm .* dRaw_F2 * dRawGrad_I1;
                                %d2SumGrad_i2   = d2SumGrad_i2 + pNorm .* dRaw_I1 * dRawGrad_F2;
                                % do FI
                                %d2SumGrad_i1   = d2SumGrad_i1 + pNorm .* dRaw_I2 * dRawGrad_F1;
                                %d2SumGrad_i2   = d2SumGrad_i2 + pNorm .* dRaw_F1 * dRawGrad_I2;
                                % do FF
                                %d2SumGrad_i1   = d2SumGrad_i1 + pNorm .* dRaw_F2 * dRawGrad_F1;
                                %d2SumGrad_i2   = d2SumGrad_i2 + pNorm .* dRaw_F1 * dRawGrad_F2;
                                
                                % do all of the gradients together
                                d2SumGrad_i1    = d2SumGrad_i1 + pNorm .* dRaw_2 * dRawGrad_1;
                                d2SumGrad_i2    = d2SumGrad_i2 + pNorm .* dRaw_1 * dRawGrad_2;
                            end
                            
                            % only do the following if we need time gradients
                            % of d2
                            if dod2TimeGrad
                                
                                % calculate gradients of mean squared doses, put
                                % them directly in d2Grad
                                % now only look at gradients wrt time
                                
                                % do II
                                %d2SumTimeGrad(temptIx_Vec)    = d2SumTimeGrad(temptIx_Vec) + pNormGrad .* d2SumRaw_II;
                                % do IF
                                %d2SumTimeGrad(temptIx_Vec)    = d2SumTimeGrad(temptIx_Vec) + pNormGrad .* d2SumRaw_IF;
                                % do FI
                                %d2SumTimeGrad(temptIx_Vec)    = d2SumTimeGrad(temptIx_Vec) + pNormGrad .* d2SumRaw_FI;
                                % do FF
                                %d2SumTimeGrad(temptIx_Vec)    = d2SumTimeGrad(temptIx_Vec) + pNormGrad .* d2SumRaw_FF;
                                
                                % combine into a single multiplication
                                d2SumTimeGrad(temptIx_Vec)    = d2SumTimeGrad(temptIx_Vec) + pNormGrad .* d2SumRaw;
                            end
                        end
                        
                        %{
                        % extract bixel weights corresponding to
                        % initial and final phases of the two beam angles
                        w_I1 = apertureInfo.arcI.bixelWeights{phase_I1}(currBixelIx1);
                        w_F1 = apertureInfo.arcF.bixelWeights{phase_F1}(currBixelIx1);
                        w_I2 = apertureInfo.arcI.bixelWeights{phase_I2}(currBixelIx2);
                        w_F2 = apertureInfo.arcF.bixelWeights{phase_F2}(currBixelIx2);
                        
                        % extract bixel jacobians corresponding to
                        % initial and final phases of the two beam angles
                        j_I1 = apertureInfo.arcI.bixelJApVec{phase_I1}(:,currBixelIx1);
                        j_F1 = apertureInfo.arcF.bixelJApVec{phase_F1}(:,currBixelIx1);
                        j_I2 = apertureInfo.arcI.bixelJApVec{phase_I2}(:,currBixelIx2);
                        j_F2 = apertureInfo.arcF.bixelJApVec{phase_F2}(:,currBixelIx2);
                        
                        % calculate gradients of mean squared doses
                        d2Grad_II = dij.scaleFactor.^2 .* ( ((dij.physicalDose{phase_I1}(dij.targetVox,currBixelIx1) * w_I1) .* (dij.physicalDose{phase_I2}(dij.targetVox,currBixelIx2) * w_I2)) * pNormGrad + ...
                            pNorm .* (dij.physicalDose{phase_I1}(dij.targetVox,currBixelIx1) * j_I1') .* repmat((dij.physicalDose{phase_I2}(dij.targetVox,currBixelIx2) * w_I2),1,size(apertureInfo.apertureVector,1)) + ...
                            pNorm .* repmat((dij.physicalDose{phase_I1}(dij.targetVox,currBixelIx1) * w_I1),1,size(apertureInfo.apertureVector,1)) .* (dij.physicalDose{phase_I2}(dij.targetVox,currBixelIx2) * j_I2') );
                        
                        d2Grad_IF = dij.scaleFactor.^2 .* ( ((dij.physicalDose{phase_I1}(dij.targetVox,currBixelIx1) * w_I1) .* (dij.physicalDose{phase_F2}(dij.targetVox,currBixelIx2) * w_F2)) * pNormGrad + ...
                            pNorm .* (dij.physicalDose{phase_I1}(dij.targetVox,currBixelIx1) * j_I1') .* repmat((dij.physicalDose{phase_F2}(dij.targetVox,currBixelIx2) * w_F2),1,size(apertureInfo.apertureVector,1)) + ...
                            pNorm .* repmat((dij.physicalDose{phase_I1}(dij.targetVox,currBixelIx1) * w_I1),1,size(apertureInfo.apertureVector,1)) .* (dij.physicalDose{phase_F2}(dij.targetVox,currBixelIx2) * j_F2') );
                        
                        d2Grad_FI = dij.scaleFactor.^2 .* ( ((dij.physicalDose{phase_F1}(dij.targetVox,currBixelIx1) * w_F1) .* (dij.physicalDose{phase_I2}(dij.targetVox,currBixelIx2) * w_I2)) * pNormGrad + ...
                            pNorm .* (dij.physicalDose{phase_F1}(dij.targetVox,currBixelIx1) * j_F1') .* repmat((dij.physicalDose{phase_I2}(dij.targetVox,currBixelIx2) * w_I2),1,size(apertureInfo.apertureVector,1)) + ...
                            pNorm .* repmat((dij.physicalDose{phase_F1}(dij.targetVox,currBixelIx1) * w_F1),1,size(apertureInfo.apertureVector,1)) .* (dij.physicalDose{phase_I2}(dij.targetVox,currBixelIx2) * j_I2') );
                        
                        d2Grad_FF = dij.scaleFactor.^2 .* ( ((dij.physicalDose{phase_F1}(dij.targetVox,currBixelIx1) * w_F1) .* (dij.physicalDose{phase_F2}(dij.targetVox,currBixelIx2) * w_F2)) * pNormGrad + ...
                            pNorm .* (dij.physicalDose{phase_F1}(dij.targetVox,currBixelIx1) * j_F1') .* repmat((dij.physicalDose{phase_F2}(dij.targetVox,currBixelIx2) * w_F2),1,size(apertureInfo.apertureVector,1)) + ...
                            pNorm .* repmat((dij.physicalDose{phase_F1}(dij.targetVox,currBixelIx1) * w_F1),1,size(apertureInfo.apertureVector,1)) .* (dij.physicalDose{phase_F2}(dij.targetVox,currBixelIx2) * j_F2') );
                        %}
                    end
                end
            end
        end
        
        % dump gradient for i2
        d2SumGrad(currVarIx2(d2KeepVar2)) = d2SumGrad(currVarIx2(d2KeepVar2)) + d2SumGrad_i2;
    end
    
    % dump gradient for i1
    d2SumGrad(currVarIx1(d2KeepVar1)) = d2SumGrad(currVarIx1(d2KeepVar1)) + d2SumGrad_i1;
end

% dump gradient for time
d2SumGrad(:,tIxBig_Vec) = d2SumGrad(:,tIxBig_Vec) + d2SumTimeGrad;

% calculate sum of variances
dVarSum = d2Sum-d'*d;

% calculate gradient of sum of variance
dVarSumGrad = d2SumGrad-2.*d'*dGrad;

end