function [dVarSum,dVarSumGrad] = matRad_doseVarianceSum(apertureInfo,dij)

%% setup

% get number of phases, subphases
numPhases       = apertureInfo.numPhases;
numSubPhases    = apertureInfo.motionModel.indices.nSubPhases;

% allocate raw doses
dRawCell_sumI   = cell(numel(apertureInfo.beam),1);
dRawCell_sumF   = cell(numel(apertureInfo.beam),1);

% allocate gradient of raw doses
dRawGradCell_sumI   = cell(numel(apertureInfo.beam),1);
dRawGradCell_sumF   = cell(numel(apertureInfo.beam),1);

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
    dRawCell_sumI{i} = zeros(numel(dij.targetVox),numSubPhases);
    dRawCell_sumF{i} = zeros(numel(dij.targetVox),numSubPhases);
    
    % allocate gradient of raw doses
    dRawGradCellTemp_sumI = zeros(numel(dij.targetVox),numel(d2KeepVar),numPhases);
    dRawGradCellTemp_sumF = zeros(numel(dij.targetVox),numel(d2KeepVar),numPhases);
    
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
                dRawCell_sumI{i}(:,subPhase_F) = dRawCell_sumI{i}(:,subPhase_F)+dRawTemp.*factor_I;
            end
            
            for subPhase_I = find(apertureInfo.motionModel.indices.subPhase2PosPhase == phase_I)'
                
                % determine factors
                factor_F = apertureInfo.probF_kL{i}(subPhase_I,phase_F)./apertureInfo.probI_IJ{i}(phase_I,phase_F);
                
                % store sums of doses (over initial and final phases)
                dRawCell_sumF{i}(:,subPhase_I) = dRawCell_sumF{i}(:,subPhase_I)+dRawTemp.*factor_F;
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
            dRawGradCellTemp_sumI(:,:,phase_F) = dRawGradCellTemp_sumI(:,:,phase_F)+dRawGradCellTemp.*factor_I;
            dRawGradCellTemp_sumF(:,:,phase_I) = dRawGradCellTemp_sumF(:,:,phase_I)+dRawGradCellTemp.*factor_F;
            
            % determine if we're doing time gradients
            deleteInd               = abs(factorGrad_F) < eps | isnan(factorGrad_F);
            factorGrad_F(deleteInd) = [];
            temptimeInd             = timeInd;
            temptimeInd(deleteInd)  = [];
            doTimeGrad              = ~all(deleteInd);
            
            if doTimeGrad
                
                [~,keepInd] = ismember(temptimeInd,d2KeepVar);
                
                dRawGradCellTemp_sumF(:,keepInd,phase_I) = dRawGradCellTemp_sumF(:,keepInd,phase_I)+dRawTemp.*factorGrad_F';
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
    
    % permute gradients, concatenate
    dRawGradCellTemp_sumI = permute(dRawGradCellTemp_sumI,[3 1 2]);
    dRawGradCellTemp_sumF = permute(dRawGradCellTemp_sumF,[1 3 2]);
    
    dRawGradCell_sumI{i} = reshape(dRawGradCellTemp_sumI,[],numel(d2KeepVar));
    dRawGradCell_sumF{i} = reshape(dRawGradCellTemp_sumF,[],numel(d2KeepVar));
    
    % dump gradients for i
    dGrad(:,currVarIx(d2KeepVar))   = dGrad(:,currVarIx(d2KeepVar))+dGrad_i;
    d2SumGrad(currVarIx(d2KeepVar)) = d2SumGrad(currVarIx(d2KeepVar))+d2SumGrad_i;
end

% first loop over all beams
for i1 = 1:numel(apertureInfo.beam)
    
    %% housekeeping for beam 1
    
    % find relevant variables
    currVarIx1 = apertureInfo.beam(i1).local2GlobalVar;
    
    % find kept variables
    d2KeepVar1 = apertureInfo.beam(i1).d2KeepVar;
    
    % get pre-calculated raw doses and gradients
    dRaw_i1     = dRawCell_sumI{i1};
    dRawGrad_i1 = dRawGradCell_sumI{i1};
    
    
    % second loop over all beams
    for i2 = (i1+1):numel(apertureInfo.beam)
        
        %% housekeeping for beam 2
        
        % find relevant variables
        currVarIx2 = apertureInfo.beam(i2).local2GlobalVar;
        
        % find kept variables
        d2KeepVar2 = apertureInfo.beam(i2).d2KeepVar;
        
        % get pre-calculated raw doses and gradients
        dRaw_i2_T   = dRawCell_sumF{i2}'; % fix this earlier?
        dRawGrad_i2 = dRawGradCell_sumF{i2};
        
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
        
        %% calculate the contribution to d2 from this combination of beam angles
        
        % get sum of dose squared
        pNormMat_x_dRaw_i2_T = (pNormMat*dRaw_i2_T)';
        pNormMat_x_dRaw_i2_T = pNormMat_x_dRaw_i2_T(:);
        d2Temp  = 2*dRaw_i1(:).'*pNormMat_x_dRaw_i2_T;
        %d2Temp  = trace(2*dRaw_i1*pNormMat*dRaw_i2);
        d2Sum   = d2Sum+d2Temp;
        
        % before we do the gradients, collapse the probability matrix along
        % the row and columns (for grad i1, i2 respectively)
        pNormMat_collapsedRows = accumarray([apertureInfo.motionModel.indices.subPhase2PosPhase_gridI(:) repelem((1:apertureInfo.motionModel.indices.nSubPhases)',apertureInfo.motionModel.indices.nSubPhases)],pNormMat(:));
        pNormMat_collapsedColumns = accumarray([repmat((1:apertureInfo.motionModel.indices.nSubPhases)',[apertureInfo.motionModel.indices.nSubPhases 1]) apertureInfo.motionModel.indices.subPhase2PosPhase_gridJ(:)],pNormMat(:));
        
        % also do some matrix mutliplication
        pNormMat_collapsed_x_dRaw_i1        = dRaw_i1*pNormMat_collapsedColumns;
        pNormMat_collapsed_x_dRaw_i1        = pNormMat_collapsed_x_dRaw_i1(:)';
        pNormMat_collapsed_x_dRaw_i2        = pNormMat_collapsedRows*dRaw_i2_T;
        pNormMat_collapsed_x_dRaw_i2        = pNormMat_collapsed_x_dRaw_i2(:)';
        pNormGradMat_x_dRaw_i2_T            = (pNormGradMat*dRaw_i2_T)';
        
        % now calculate gradients
        d2SumGrad_i1        = 2*pNormMat_collapsed_x_dRaw_i2*dRawGrad_i1;
        d2SumGrad_i2        = 2*pNormMat_collapsed_x_dRaw_i1*dRawGrad_i2;
        d2SumTimeGrad_i1i2  = 2*dRaw_i1(:).'*pNormGradMat_x_dRaw_i2_T(:);
        
        % dump gradient for beams 1 and 2
        d2SumGrad(currVarIx1(d2KeepVar1))   = d2SumGrad(currVarIx1(d2KeepVar1)) + d2SumGrad_i1;
        d2SumGrad(currVarIx2(d2KeepVar2))   = d2SumGrad(currVarIx2(d2KeepVar2)) + d2SumGrad_i2;
        d2SumTimeGrad                       = d2SumTimeGrad + d2SumTimeGrad_i1i2*pNormGradFactor;
    end
end

% dump gradient for time
d2SumGrad(:,tIxBig_Vec) = d2SumGrad(:,tIxBig_Vec) + d2SumTimeGrad;

% calculate sum of variances
dVarSum = d2Sum-d'*d;

% calculate gradient of sum of variance
dVarSumGrad = d2SumGrad-2.*d'*dGrad;

end