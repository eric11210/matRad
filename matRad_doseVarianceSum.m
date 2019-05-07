function [dVarSum,dVarSumGrad] = matRad_doseVarianceSum(apertureInfo,dij)

%% setup

% INSERT AND EXTRA FACTOR OF numPhases HERE

% allocate raw doses
dRaw    = cell(numel(apertureInfo.beam).*apertureInfo.numPhases.^2,1);
%dRaw_F = cell(numel(apertureInfo.beam).*apertureInfo.numPhases.^2,1);

% allocate gradient of raw doses
dRawGradCell    = cell(numel(apertureInfo.beam).*apertureInfo.numPhases.^2,1);
%dRawGradCell_F = cell(numel(apertureInfo.beam).*apertureInfo.numPhases.^2,1);

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
nnzGrad     = 0;
%nnzGrad_F  = 0;

for i = 1:numel(apertureInfo.beam)
    
    % pre-calculate raw doses
    % first find relevant bixels
    currBixelIx = apertureInfo.beam(i).bixelIndMap(~isnan(apertureInfo.beam(i).bixelIndMap));
    
    % find relevant variables for beam
    currVarIx = apertureInfo.beam(i).local2GlobalVar;
    
    % allocate gradient of dose for beam
    dGrad_i = zeros(numel(dij.targetVox),apertureInfo.beam(i).numUniqueVar);
    
    % loop over initial and final phases
    for phase_I = 1:apertureInfo.numPhases
        
        % extract dijBeamPhase for phase_I
        dijBeamPhase_I      = dij.scaleFactor .* dij.physicalDose{phase_I}(dij.targetVox,currBixelIx);
        dijBeamPhaseInv_I   = dijBeamPhase_I';
        
        for phase_F = 1:apertureInfo.numPhases
            
            % determine cell index
            cellInd = (phase_I-1).*apertureInfo.numPhases+phase_F;
            
            % extract dijBeamPhase for phase_F
            if phase_I == phase_F
                dijBeamPhase_F      = dijBeamPhase_I;
                dijBeamPhaseInv_F   = dijBeamPhaseInv_I;
            else
                dijBeamPhase_F      = dij.scaleFactor .* dij.physicalDose{phase_F}(dij.targetVox,currBixelIx);
                dijBeamPhaseInv_F   = dijBeamPhase_F';
            end
            
            % now extract bixel weights
            w_I = apertureInfo.arcI.bixelWeights{cellInd}(currBixelIx)';
            w_F = apertureInfo.arcF.bixelWeights{cellInd}(currBixelIx)';
            % calculate dose
            dRawTemp_I  = w_I * dijBeamPhaseInv_I;
            dRawTemp_F  = w_F * dijBeamPhaseInv_F;
            dRawTemp    = dRawTemp_I+dRawTemp_F;
            
            % calculate gradient of dose
            % d = d + dij.scaleFactor .*
            % dij.physicalDose{phase_I1}(:,currBixelIx1) * apertureInfo.bixelWeights{phase_I1}(currBixelIx1);
            j_I = apertureInfo.arcI.bixelJApVec{(i-1).*apertureInfo.numPhases.^2+cellInd}';
            j_F = apertureInfo.arcF.bixelJApVec{(i-1).*apertureInfo.numPhases.^2+cellInd}';
            % calculate gradient
            dRawGradCellTemp_I  = dijBeamPhase_I * j_I;
            dRawGradCellTemp_F  = dijBeamPhase_F * j_F;
            dRawGradCellTemp    = dRawGradCellTemp_I+dRawGradCellTemp_F;
            
            % accumulate sum of dose and dose gradient
            d       = d+dRawTemp';
            dGrad_i = dGrad_i+dRawGradCellTemp;
            
            % store doses and gradients
            dRaw{(i-1).*apertureInfo.numPhases.^2+cellInd}          = dRawTemp;
            dRawGradCell{(i-1).*apertureInfo.numPhases.^2+cellInd}  = dRawGradCellTemp;
            %dRaw_F{(i-1).*apertureInfo.numPhases.^2+cellInd} = dRawTemp_F;
            %dRawGradCell_F{(i-1).*apertureInfo.numPhases.^2+cellInd} = dRawGradCellTemp_F;
            
            % track number of nonzero elements
            nnzGrad     = nnzGrad+nnz(dRawGradCellTemp);
            %nnzGrad_F  = nnzGrad_F+nnz(dRawGradCellTemp_F);
            
        end
    end
    
    % dump gradient for i
    dGrad(:,currVarIx) = dGrad(:,currVarIx)+dGrad_i;
end

% preallocate full values and indices
V = zeros(nnzGrad,1);
I = zeros(nnzGrad,1);
J = zeros(nnzGrad,1);
%V_F = zeros(nnzGrad_F,1);
%I_F = zeros(nnzGrad_F,1);
%J_F = zeros(nnzGrad_F,1);
%dRawGrad_I = spalloc(apertureInfo.beam(end).gradOffset(end)+apertureInfo.beam(end).numUniqueVar-1,numel(dij.targetVox),nnzGrad_I);
%dRawGrad_F = spalloc(apertureInfo.beam(end).gradOffset(end)+apertureInfo.beam(end).numUniqueVar-1,numel(dij.targetVox),nnzGrad_F);

% initialize offset
indOffset = 0;
%indOffset_F = 0;

% loop over all beams again to put individual values and indices into full
for i = 1:numel(apertureInfo.beam)
    
    % get gradient offset
    gradOffset = apertureInfo.beam(i).gradOffset;
    
    % get number of unique variables
    numUniqueVar = apertureInfo.beam(i).numUniqueVar;
    
    % loop over phases
    for phase_I = 1:apertureInfo.numPhases
        
        for phase_F = 1:apertureInfo.numPhases
            
            % determine cell index
            cellInd = (phase_I-1).*apertureInfo.numPhases+phase_F;
            
            % extract gradients
            dRawGradTemp    = dRawGradCell{(i-1).*apertureInfo.numPhases.^2+cellInd};
            %dRawGradTemp_F = dRawGradCell_F{(i-1).*apertureInfo.numPhases.^2+cellInd};
            
            % get nonzero values and indices
            vTemp   = nonzeros(dRawGradTemp);
            %v_F    = nonzeros(dRawGradTemp_F);
            [iTemp,jTemp]   = find(dRawGradTemp);
            %[i_F,j_F]      = find(dRawGradTemp_F);
            
            % get indices for full
            ind     = indOffset+(1:numel(vTemp));
            %ind_F  = indOffset_F+(1:numel(v_F));
            
            % put values and indives into full
            V(ind)      = vTemp;
            I(ind)      = iTemp;
            J(ind)      = jTemp+gradOffset+(cellInd-1).*numUniqueVar;
            %V_F(ind_F) = v_F;
            %I_F(ind_F) = i_F;
            %J_F(ind_F) = j_F+gradOffset+(cellInd-1).*numUniqueVar;
            
            % wipe cell
            dRawGradCell{(i-1).*apertureInfo.numPhases.^2+cellInd}      = [];
            %dRawGradCell_F{(i-1).*apertureInfo.numPhases.^2+cellInd}   = [];
            
            % increase offset
            indOffset       = indOffset+numel(vTemp);
            %indOffset_F    = indOffset_F+numel(v_F);
            
        end
    end
end

% construct sparse matrix of gradients
dRawGrad = sparse(I,J,V,numel(dij.targetVox),apertureInfo.beam(end).gradOffset+apertureInfo.beam(end).numUniqueVar.*apertureInfo.numPhases.^2,nnzGrad);
%dRawGrad_F = sparse(I_F,J_F,V_F,numel(dij.targetVox),apertureInfo.beam(end).gradOffset+apertureInfo.beam(end).numUniqueVar.*apertureInfo.numPhases.^2,nnzGrad_F);

%%% DO SUMMATION OF I AND F FARTHER UP? IS THERE ANY NEED TO KEEP THEM
%%% SEPARATE AT THIS POINT?

% clear unnecessary variables
clear I_I J_I V_I I_F J_F V_F

% first loop over all beams
for i1_loop = 1:numel(apertureInfo.beam)
    
    % second loop over all beams
    for i2_loop = 1:numel(apertureInfo.beam)
        
        %% housekeeping
        
        % determine the real i1 and i2 from the loop variables
        % hint: i1 must be smaller than or equal to i2
        i1  = min([i1_loop i2_loop]);
        i2  = max([i1_loop i2_loop]);
        
        % find relevant variables for beams 1 and 2
        currVarIx1 = apertureInfo.beam(i1).local2GlobalVar;
        currVarIx2 = apertureInfo.beam(i2).local2GlobalVar;
        
        % find kept variables for beams 1 and 2
        d2KeepVar1 = apertureInfo.beam(i1).d2KeepVar;
        d2KeepVar2 = apertureInfo.beam(i2).d2KeepVar;
        %d2KeepVar1(:) = true;
        %d2KeepVar1(6:end) = false;
        
        % find variable indices for beams 1 and 2
        numUniqueVar_i1 = apertureInfo.beam(i1).numUniqueVar;
        varInd_i1 = apertureInfo.beam(i1).gradOffset+(1:apertureInfo.beam(i1).numUniqueVar);
        varInd_i1(~d2KeepVar1) = [];
        numUniqueVar_i2 = apertureInfo.beam(i2).numUniqueVar;
        varInd_i2 = apertureInfo.beam(i2).gradOffset+(1:apertureInfo.beam(i2).numUniqueVar);
        varInd_i2(~d2KeepVar2) = [];
        
        % allocate gradient of mean of squared dose for beams 1 and 2
        d2SumGrad_i1 = zeros(1,nnz(d2KeepVar1));
        d2SumGrad_i2 = zeros(1,nnz(d2KeepVar2));
        %d2SumGrad_i1 = zeros(1,apertureInfo.beam(i1).numUniqueVar);
        
        %% determine probabilities
        % 1) probability to transfer from the end of the earlier to the
        % beginning of the later
        % 2) probability to arrive at the later
        
        % determine pNormGradFactors
        pNormGradFactor1 = sum(apertureInfo.propVMAT.jacobT(:,(i1+1):(i2-1)),2);
        pNormGradFactor2 = sum(apertureInfo.propVMAT.jacobT(:,1:(i1-1)),2);
        
        % calculate the transition and arrival times
        transT12    = sum([apertureInfo.beam((i1+1):(i2-1)).time]);
        T2          = sum([apertureInfo.beam(1:(i2-1)).time]);
        
        % calculate the probabilities
        [Pij_transT12,Pij_transT12_dot,Pi_T2,Pi_T2_dot] = matRad_transAndTProb(transT12,T2,apertureInfo.motionModel);
        
        % first loop over initial phases
        for phase_I1 = 1:apertureInfo.numPhases
            
            % first loop over final phases
            for phase_F1 = 1:apertureInfo.numPhases
                
                % determine cell index for beam 1
                cellInd_1 = (phase_I1-1).*apertureInfo.numPhases+phase_F1;
                
                % get pre-calculated raw doses
                dRaw_1      = dRaw{(i1-1).*apertureInfo.numPhases+cellInd_1};
                %dRaw_F1    = dRaw_F{(i1-1).*apertureInfo.numPhases+cellInd_1};
                
                % get pre-calculated raw gradients
                % only consider gradients wrt weights and times for d2
                %varInd     = (phase_I1-1).*numUniqueVar_i1+varInd_i1;
                varInd      = (cellInd_1-1).*numUniqueVar_i1+varInd_i1;
                dRawGrad_1  = dRawGrad(:,varInd);
                %dRawGrad_F1 = dRawGrad_F(:,varInd);
                %dRawGrad_I1 = dRawGradCell_I{(i1-1).*apertureInfo.numPhases+phase_I1}(:,d2KeepVar1);
                %dRawGrad_F1 = dRawGradCell_F{(i1-1).*apertureInfo.numPhases+phase_F1}(:,d2KeepVar1);
                
                % sum I1 and F1 doses and gradients
                %dRaw_1      = dRaw_I1+dRaw_F1;
                %dRawGrad_1  = dRawGrad_I1+dRawGrad_F1;
                
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
                    
                    % second loop over final phases
                    for phase_F2 = 1:apertureInfo.numPhases
                        
                        % determine cell index for beam 2
                        cellInd_2 = (phase_I2-1).*apertureInfo.numPhases+phase_F2;
                        
                        %% determine mean of squared dose for this combination
                        
                        % check if we are doing either d2 or dod2TimeGrad stuff
                        if dod2 || dod2TimeGrad
                            
                            % get pre-calculated raw doses
                            dRaw_2      = dRaw{(i2-1).*apertureInfo.numPhases+cellInd_2};
                            %dRaw_F2    = dRaw_F{(i2-1).*apertureInfo.numPhases+cellInd_2};
                            
                            % sum I2 and F2 doses
                            %dRaw_2 = dRaw_I2+dRaw_F2;
                            
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
                                %varInd     = (phase_I2-1).*numUniqueVar_i2+varInd_i2;
                                varInd      = (cellInd_2-1).*numUniqueVar_i2+varInd_i2;
                                dRawGrad_2  = dRawGrad(:,varInd);
                                %dRawGrad_F2 = dRawGrad_F(:,varInd);
                                %dRawGrad_I2 = dRawGradCell_I{(i2-1).*apertureInfo.numPhases+phase_I2}(:,d2KeepVar2);
                                %dRawGrad_F2 = dRawGradCell_F{(i2-1).*apertureInfo.numPhases+phase_F2}(:,d2KeepVar2);
                                
                                % sum I2 and F2 gradients
                                %dRawGrad_2  = dRawGrad_I2+dRawGrad_F2;
                                
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
        
        % dump gradient for beams 1 and 2
        d2SumGrad(currVarIx1(d2KeepVar1)) = d2SumGrad(currVarIx1(d2KeepVar1)) + d2SumGrad_i1;
        d2SumGrad(currVarIx2(d2KeepVar2)) = d2SumGrad(currVarIx2(d2KeepVar2)) + d2SumGrad_i2;
    end
    
    1;
end

% dump gradient for time
d2SumGrad(:,tIxBig_Vec) = d2SumGrad(:,tIxBig_Vec) + d2SumTimeGrad;

% calculate sum of variances
dVarSum = d2Sum-d'*d;

% calculate gradient of sum of variance
dVarSumGrad = d2SumGrad-2.*d'*dGrad;

end