function [dVarSum,dVarSumGrad] = matRad_doseVarianceSum_SLOW(apertureInfo,dij)

%% setup

% get number of phases, subphases
numPhases       = apertureInfo.numPhases;
numSubPhases    = apertureInfo.motionModel.indices.nSubPhases;

% allocate raw doses
dRawCell_sumI   = cell(numel(apertureInfo.beam).*numSubPhases,1);
dRawCell_sumF   = cell(numel(apertureInfo.beam).*numSubPhases,1);

% allocate gradient of raw doses
dRawGradCell_sumI   = cell(numel(apertureInfo.beam).*numPhases,1);
dRawGradCell_sumF   = cell(numel(apertureInfo.beam).*numPhases,1);

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
tIxBig_Vec  = (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*numPhases+(1:apertureInfo.totalNumOfShapes);
tIx_Vec     = 1:apertureInfo.totalNumOfShapes;

%% precalculate doses and gradients
% loop over all beams

% track nnz
%nnzGrad     = 0;

for i = 1:numel(apertureInfo.beam)
    
    % pre-calculate raw doses
    % first find relevant bixels
    currBixelIx = apertureInfo.beam(i).bixelIndMap(~isnan(apertureInfo.beam(i).bixelIndMap));
    
    % find relevant variables for beam
    currVarIx = apertureInfo.beam(i).local2GlobalVar;
    
    % find kept variables
    d2KeepVar = apertureInfo.beam(i).d2KeepVar;
    
    % allocate gradient of dose and square of dose for beam
    dGrad_i     = zeros(numel(dij.targetVox),apertureInfo.beam(i).numKeepVar);
    d2SumGrad_i = zeros(1,apertureInfo.beam(i).numKeepVar);
    
    % allocate raw doses
    dRawCell_sumI((i-1).*numSubPhases+(1:numSubPhases)) = {zeros(1,numel(dij.targetVox))};
    dRawCell_sumF((i-1).*numSubPhases+(1:numSubPhases)) = {zeros(1,numel(dij.targetVox))};
    
    % allocate gradient of raw doses
    dRawGradCell_sumI((i-1).*numPhases+(1:numPhases)) = {zeros(numel(dij.targetVox),numel(d2KeepVar))};
    dRawGradCell_sumF((i-1).*numPhases+(1:numPhases)) = {zeros(numel(dij.targetVox),numel(d2KeepVar))};
    
    % calculate probability normalization and gradient matrices
    pNormMat        = 1./apertureInfo.probI_IJ{i};
    pNormGradMat    = -apertureInfo.probIGrad_IJ{i}./(apertureInfo.probI_IJ{i}.^2);
    
    % time indices
    timeInd = tIx_Vec+apertureInfo.beam(i).numUniqueVar-numel(apertureInfo.beam);
    
    % loop over initial and final phases
    for phase_I = 1:numPhases
        
        % extract dijBeamPhase for phase_I
        dijBeamPhase_I      = dij.scaleFactor .* dij.physicalDose{phase_I}(dij.targetVox,currBixelIx);
        
        for phase_F = 1:numPhases
            
            % determine cell index
            cellInd = (phase_I-1).*numPhases+phase_F;
            
            % extract dijBeamPhase for phase_F
            if phase_I == phase_F
                dijBeamPhase_F      = dijBeamPhase_I;
            else
                dijBeamPhase_F      = dij.scaleFactor .* dij.physicalDose{phase_F}(dij.targetVox,currBixelIx);
            end
            
            % now extract bixel weights
            w_I = apertureInfo.arcI.bixelWeights{cellInd}(currBixelIx);
            w_F = apertureInfo.arcF.bixelWeights{cellInd}(currBixelIx);
            % calculate dose
            dRawTemp = dijBeamPhase_I * w_I+dijBeamPhase_F * w_F;
            
            % calculate gradient of dose
            % d = d + dij.scaleFactor .*
            % dij.physicalDose{phase_I1}(:,currBixelIx1) * apertureInfo.bixelWeights{phase_I1}(currBixelIx1);
            j_I = apertureInfo.arcI.bixelJApVec{(i-1).*numPhases.^2+cellInd}';
            j_F = apertureInfo.arcF.bixelJApVec{(i-1).*numPhases.^2+cellInd}';
            
            j_I = j_I(:,d2KeepVar);
            j_F = j_F(:,d2KeepVar);
            % calculate gradient
            dRawGradCellTemp    = dijBeamPhase_I * j_I+dijBeamPhase_F * j_F;
            
            % accumulate sum of dose and dose gradient
            d       = d+dRawTemp;
            dGrad_i = dGrad_i+dRawGradCellTemp;
            
            % track number of nonzero elements
            %nnzGrad     = nnzGrad+nnz(dRawGradCellTemp);
            
            % store sums of doses and gradients (over initial and final
            % phases)
            %dRawCell_sumI{(i-1).*numPhases+phase_F} = dRawCell_sumI{(i-1).*numPhases+phase_F}+dRawTemp';
            %dRawCell_sumF{(i-1).*numPhases+phase_I} = dRawCell_sumF{(i-1).*numPhases+phase_I}+dRawTemp';
            
            for subPhase_F = find(apertureInfo.motionModel.indices.subPhase2PosPhase == phase_F)'
                
                % determine factors
                factor_I = apertureInfo.probI_Ij{i}(phase_I,subPhase_F)./apertureInfo.probI_IJ{i}(phase_I,phase_F);
                
                % store sums of doses (over initial and final phases)
                dRawCell_sumI{(i-1).*numSubPhases+subPhase_F} = dRawCell_sumI{(i-1).*numSubPhases+subPhase_F}+dRawTemp'.*factor_I;
            end
            
            for subPhase_I = find(apertureInfo.motionModel.indices.subPhase2PosPhase == phase_I)'
                
                % determine factors
                factor_F = apertureInfo.probF_kL{i}(subPhase_I,phase_F)./apertureInfo.probI_IJ{i}(phase_I,phase_F);
                
                % store sums of doses (over initial and final phases)
                dRawCell_sumF{(i-1).*numSubPhases+subPhase_I} = dRawCell_sumF{(i-1).*numSubPhases+subPhase_I}+dRawTemp'.*factor_F;
            end
            
            % determine factors and gradients for sums
            factor_I        = 1./apertureInfo.motionModel.indices.nSubPhasePerPosPhase(phase_F);
            %factorGrad_I    = 0;
            factor_F        = apertureInfo.probF_KL{i}(phase_I,phase_F)./(apertureInfo.probI_IJ{i}(phase_I,phase_F).*apertureInfo.motionModel.indices.nSubPhasePerPosPhase(phase_I));
            factorGrad_F    = squeeze((1./apertureInfo.motionModel.indices.nSubPhasePerPosPhase(phase_I)).*(apertureInfo.probFGrad_KL{i}(phase_I,phase_F,:)./apertureInfo.probI_IJ{i}(phase_I,phase_F) - ...
                apertureInfo.probF_KL{i}(phase_I,phase_F).*apertureInfo.probIGrad_IJ{i}(phase_I,phase_F,:)./(apertureInfo.probI_IJ{i}(phase_I,phase_F).^2)));
            
            % store sums of doses (over initial and final phases)
            %dRawCell_sumI{(i-1).*numPhases+phase_F} = dRawCell_sumI{(i-1).*numPhases+phase_F}+dRawTemp'.*factor_I;
            %dRawCell_sumF{(i-1).*numPhases+phase_I} = dRawCell_sumF{(i-1).*numPhases+phase_I}+dRawTemp'.*factor_F;
            
            % insert gradients (averaged over all sub-phases in a phase
            dRawGradCell_sumI{(i-1).*numPhases+phase_F} = dRawGradCell_sumI{(i-1).*numPhases+phase_F}+dRawGradCellTemp.*factor_I;
            dRawGradCell_sumF{(i-1).*numPhases+phase_I} = dRawGradCell_sumF{(i-1).*numPhases+phase_I}+dRawGradCellTemp.*factor_F;
            
            % determine if we're doing time gradients
            deleteInd               = abs(factorGrad_F) < eps | isnan(factorGrad_F);
            factorGrad_F(deleteInd) = [];
            temptimeInd             = timeInd;
            temptimeInd(deleteInd)  = [];
            doTimeGrad              = ~all(deleteInd);
            
            if doTimeGrad
                
                [~,keepInd] = ismember(temptimeInd,d2KeepVar);
                
                dRawGradCell_sumF{(i-1).*numPhases+phase_I}(:,keepInd) = dRawGradCell_sumF{(i-1).*numPhases+phase_I}(:,keepInd)+dRawTemp.*factorGrad_F';
            end
            
            
            %% do d2 for i1 == i2 now
            
            % determine probability normalization and gradient for d2
            pNorm       = pNormMat(phase_I,phase_F);
            pNormGrad   = squeeze(pNormGradMat(phase_I,phase_F,:));
            
            % determine if we're doing d2 for this combo
            dod2 = abs(pNorm) > eps && pNorm < inf;
            
            % determine if we're doing time gradients of d2 for
            % this combo
            deleteInd                   = abs(pNormGrad) < eps | isnan(pNormGrad);
            pNormGrad(deleteInd)        = [];
            temptIx_Vec                 = tIx_Vec;
            temptIx_Vec(deleteInd)      = [];
            dod2TimeGrad                = ~all(deleteInd);
            
            % check if we are doing either d2 or dod2TimeGrad stuff
            if dod2 || dod2TimeGrad
                
                if dod2
                    
                    % combine into a single multiplication
                    d2SumRaw = dRawTemp'*dRawTemp;
                    
                    % sum terms
                    d2Sum = d2Sum + pNorm .* d2SumRaw;
                    
                    % do matrix multiplication once to save time
                    addedGrad       = 2 .* pNorm .* dRawTemp' * dRawGradCellTemp;
                    d2SumGrad_i    = d2SumGrad_i + addedGrad;
                end
                
                % only do the following if we need time gradients
                % of d2
                if dod2TimeGrad
                    
                    % combine into a single multiplication
                    d2SumTimeGrad(temptIx_Vec)    = d2SumTimeGrad(temptIx_Vec) + pNormGrad' .* d2SumRaw;
                end
            end
        end
    end
    
    % dump gradients for i
    dGrad(:,currVarIx(d2KeepVar))   = dGrad(:,currVarIx(d2KeepVar))+dGrad_i;
    d2SumGrad(currVarIx(d2KeepVar)) = d2SumGrad(currVarIx(d2KeepVar))+d2SumGrad_i;
end

%{
% preallocate full values and indices
V = zeros(nnzGrad,1);
I = zeros(nnzGrad,1);
J = zeros(nnzGrad,1);

% initialize offset
indOffset = 0;

% extract gantry times
gantryTimes = zeros(numel(apertureInfo.beam),1);

% loop over all beams again to put individual values and indices into full
for i = 1:numel(apertureInfo.beam)
    
    % get gradient offset
    gradOffset = apertureInfo.beam(i).gradOffset;
    
    % get number of kept variables
    numKeepVar = apertureInfo.beam(i).numKeepVar;
    
    % extract times
    gantryTimes(i) = apertureInfo.beam(i).time;
    
    % loop over phases
    for phase_I = 1:numPhases
        
        for phase_F = 1:numPhases
            
            % determine cell index
            cellInd = (phase_I-1).*numPhases+phase_F;
            
            % extract gradients
            dRawGradTemp    = dRawGradCell{(i-1).*numPhases.^2+cellInd};
            
            % get nonzero values and indices
            vTemp   = nonzeros(dRawGradTemp);
            [iTemp,jTemp]   = find(dRawGradTemp);
            
            % get indices for full
            ind     = indOffset+(1:numel(vTemp));
            
            % put values and indives into full
            V(ind)      = vTemp;
            I(ind)      = iTemp;
            J(ind)      = jTemp+gradOffset+(cellInd-1).*numKeepVar;
            
            % increase offset
            indOffset       = indOffset+numel(vTemp);
        end
    end
end

% construct sparse matrix of gradients
dRawGrad = sparse(I,J,V,numel(dij.targetVox),apertureInfo.beam(end).gradOffset+apertureInfo.beam(end).numKeepVar.*numPhases.^2);

% clear unnecessary variables
clear I J V
%}
% first loop over all beams
for i1 = 1:numel(apertureInfo.beam)
    
    % second loop over all beams
    for i2 = (i1+1):numel(apertureInfo.beam)
        
        %% housekeeping
        
        % find relevant variables for beams 1 and 2
        currVarIx1 = apertureInfo.beam(i1).local2GlobalVar;
        currVarIx2 = apertureInfo.beam(i2).local2GlobalVar;
        
        % find kept variables for beams 1 and 2
        d2KeepVar1 = apertureInfo.beam(i1).d2KeepVar;
        d2KeepVar2 = apertureInfo.beam(i2).d2KeepVar;
        
        % allocate gradient of mean of squared dose for beams 1 and 2
        d2SumGrad_i1 = zeros(1,nnz(d2KeepVar1));
        d2SumGrad_i2 = zeros(1,nnz(d2KeepVar2));
        
        %% determine probabilities
        
        % determine probability normalization matrix
        % if we are looking at different beams, we need to normalize by
        % multiplying by (the probability to transition from F1 to I2)
        % divided (by the probability to arrive at I2)
        
        % determine pNormGradFactor
        pNormGradFactor = sum(apertureInfo.propVMAT.jacobT(:,(i1+1):(i2-1)),2)';
        
        % calculate the transition and arrival times
        transT12    = sum([apertureInfo.beam((i1+1):(i2-1)).time]);
        
        % calculate the probabilities
        [pNormMat,pNormGradMat,~,~] = matRad_transAndTProb(transT12,0,apertureInfo.motionModel);
        
        % determine thresholds
        pNormThres      = eps;%max(pNormMat(:))*1e-2;
        pNormGradThres  = eps;%max(pNormGradMat(:))*1e-2;
        
        %% loop over phases F1 and I2 only; we have already performed a sum over I1 and F2
        % if we are looking at different beams, then we need to
        % loop over initial and final phases for the second
        % beam
        
        % first loop over final phases
        for phase_F1 = 1:numPhases
            
            % find subphases in phase_F1
            subPhases_F1 = find(apertureInfo.motionModel.indices.subPhase2PosPhase == phase_F1)';
            
            % get pre-calculated raw gradients
            dRawGrad_1 = dRawGradCell_sumI{(i1-1).*numPhases+phase_F1};
            
            % second loop over initial phases
            for phase_I2 = 1:numPhases
                
                % find subphases in phase_I2
                subPhases_I2 = find(apertureInfo.motionModel.indices.subPhase2PosPhase == phase_I2)';
                
                % get pre-calculated raw gradients
                dRawGrad_2 = dRawGradCell_sumF{(i2-1).*numPhases+phase_I2};
                
                % initialize sums of doses
                dSum_1 = zeros(1,numel(dij.targetVox));
                dSum_2 = zeros(1,numel(dij.targetVox));
                
                % loop over subphases in phase_F1
                for subPhase_F1 = subPhases_F1
                    
                    % get pre-calculated raw doses
                    dRaw_1 = dRawCell_sumI{(i1-1).*numSubPhases+subPhase_F1};
                    
                    % loop over subphases in phase_I2
                    for subPhase_I2 = subPhases_I2
                        
                        % determine probability normalization
                        pNorm = pNormMat(subPhase_F1,subPhase_I2);
                        
                        % determine gradient of probability normalization
                        pNormGrad = pNormGradFactor.*pNormGradMat(subPhase_F1,subPhase_I2);
                        
                        % determine if we're doing d2 for this combo
                        dod2 = abs(pNorm) > pNormThres && pNorm < inf;
                        
                        % determine if we're doing time gradients of d2 for
                        % this combo
                        %keepInd         = abs(pNormGrad) > pNormGradThres & ~isnan(pNormGrad);
                        dod2TimeGrad    = any(abs(pNormGrad) > pNormGradThres & ~isnan(pNormGrad));
                        
                        % only do the following if we need d2 stuff,
                        % i.e. if pNorm is > eps
                        if dod2 || dod2TimeGrad
                            
                            % get pre-calculated raw doses
                            dRaw_2      = dRawCell_sumF{(i2-1).*numSubPhases+subPhase_I2};
                            
                            % combine into a single multiplication
                            d2SumRaw = dRaw_1 * dRaw_2';
                            
                            % determine if we're doing d2 for this combo
                            if dod2
                                
                                % sum terms
                                d2Sum = d2Sum + 2 .* pNorm .* d2SumRaw;
                                
                                % keep track of the sum of for the purposes
                                % of gradient
                                dSum_1 = dSum_1 + 2 .* pNorm .* dRaw_1;
                                dSum_2 = dSum_2 + 2 .* pNorm .* dRaw_2;
                                
                                %d2SumGrad_i1    = d2SumGrad_i1 + 2.*pNorm .* dRaw_2 * dRawGrad_1;
                                %d2SumGrad_i2    = d2SumGrad_i2 + 2.*pNorm .* dRaw_1 * dRawGrad_2;
                            end
                            
                            % only do the following if we need time gradients
                            % of d2
                            if dod2TimeGrad
                                
                                % calculate gradients of mean squared doses, put
                                % them directly in d2Grad
                                % now only look at gradients wrt time
                                
                                % cut indices
                                %pNormGrad_cut   = pNormGrad(keepInd);
                                %tIx_Vec_cut     = tIx_Vec(keepInd);
                                
                                % combine into a single multiplication
                                d2SumTimeGrad    = d2SumTimeGrad + 2.* pNormGrad .* d2SumRaw;
                            end
                        end
                    end
                end
                
                % do all of the gradients together
                d2SumGrad_i1    = d2SumGrad_i1 + dSum_2 * dRawGrad_1;
                d2SumGrad_i2    = d2SumGrad_i2 + dSum_1 * dRawGrad_2;
            end
        end
        
        % dump gradient for beams 1 and 2
        d2SumGrad(currVarIx1(d2KeepVar1)) = d2SumGrad(currVarIx1(d2KeepVar1)) + d2SumGrad_i1;
        d2SumGrad(currVarIx2(d2KeepVar2)) = d2SumGrad(currVarIx2(d2KeepVar2)) + d2SumGrad_i2;
    end
end

% dump gradient for time
d2SumGrad(:,tIxBig_Vec) = d2SumGrad(:,tIxBig_Vec) + d2SumTimeGrad;

% calculate sum of variances
dVarSum = d2Sum-d'*d;

% calculate gradient of sum of variance
dVarSumGrad = d2SumGrad-2.*d'*dGrad;

end