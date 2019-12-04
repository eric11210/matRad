function dVar = matRad_doseVariance(apertureInfo,dij)

%% setup

% get number of phases, subphases
numPhases       = apertureInfo.numPhases;
numSubPhases    = apertureInfo.motionModel.indices.nSubPhases;

% allocate raw doses
dRawCell_sumI   = cell(numel(apertureInfo.beam).*numSubPhases,1);
dRawCell_sumF   = cell(numel(apertureInfo.beam).*numSubPhases,1);

% allocate mean of dose
d = zeros(numel(dij.targetVox),1);

% allocate mean of squared dose
d2 = zeros(numel(dij.targetVox),1);

%% precalculate doses
% loop over all beams
for i = 1:numel(apertureInfo.beam)
    
    % pre-calculate raw doses
    % get bixel and fixel indices for last/next dose beams
    lastBixelIx = apertureInfo.beam(i).lastBixelIndMap(~isnan(apertureInfo.beam(i).lastBixelIndMap));
    nextBixelIx = apertureInfo.beam(i).nextBixelIndMap(~isnan(apertureInfo.beam(i).nextBixelIndMap));
    lastFixelIx = apertureInfo.beam(i).lastFixelIndMap(~isnan(apertureInfo.beam(i).lastFixelIndMap));
    nextFixelIx = apertureInfo.beam(i).nextFixelIndMap(~isnan(apertureInfo.beam(i).nextFixelIndMap));
    
    % allocate raw doses
    dRawCell_sumI((i-1).*numSubPhases+(1:numSubPhases)) = {zeros(numel(dij.targetVox),1)};
    dRawCell_sumF((i-1).*numSubPhases+(1:numSubPhases)) = {zeros(numel(dij.targetVox),1)};
    
    % calculate probability normalization matrix
    pNormMat = 1./apertureInfo.probI_IJ{i};
    
    % loop over initial and final phases
    for phase_I = 1:numPhases
        
        % extract dij for last dose beam and phase_I
        dij_lastDose_phaseI = dij.scaleFactor .* dij.physicalDose{phase_I}(dij.targetVox,lastBixelIx);
        
        % extract dij for next dose beam and phase_I
        if all(lastBixelIx == nextBixelIx)
            dij_nextDose_phaseI = dij_lastDose_phaseI;
        else
            dij_nextDose_phaseI = dij.scaleFactor .* dij.physicalDose{phase_I}(dij.targetVox,nextBixelIx);
        end
        
        for phase_F = 1:numPhases
            
            % determine cell index
            cellInd = (phase_I-1).*numPhases+phase_F;
            
            % extract dij for last dose beam and phase_F
            if phase_I == phase_F
                dij_lastDose_phaseF      = dij_lastDose_phaseI;
            else
                dij_lastDose_phaseF      = dij.scaleFactor .* dij.physicalDose{phase_F}(dij.targetVox,lastBixelIx);
            end
            
            % extract dij for next dose beam and phase_F
            if all(lastBixelIx == nextBixelIx)
                dij_nextDose_phaseF = dij_lastDose_phaseF;
            else
                dij_nextDose_phaseF = dij.scaleFactor .* dij.physicalDose{phase_F}(dij.targetVox,nextBixelIx);
            end
            
            % now extract fixel weights for last/next dose beam and
            % arc_I/arc_F
            w_lastDose_arcI = apertureInfo.arcI.lastDose.fixelWeights{cellInd}(lastFixelIx);
            w_lastDose_arcF = apertureInfo.arcF.lastDose.fixelWeights{cellInd}(lastFixelIx);
            w_nextDose_arcI = apertureInfo.arcI.nextDose.fixelWeights{cellInd}(nextFixelIx);
            w_nextDose_arcF = apertureInfo.arcF.nextDose.fixelWeights{cellInd}(nextFixelIx);
            % calculate dose
            dRawTemp = dij_lastDose_phaseI*w_lastDose_arcI + dij_lastDose_phaseF*w_lastDose_arcF + dij_nextDose_phaseI*w_nextDose_arcI + dij_nextDose_phaseF*w_nextDose_arcF;
            
            % accumulate sum of dose
            d = d+dRawTemp;
            
            % determine factors for sums
            %factor_F = apertureInfo.probF_KL{i}(phase_I,phase_F)./(apertureInfo.probI_IJ{i}(phase_I,phase_F).*apertureInfo.motionModel.indices.nSubPhasePerPosPhase(phase_I));
            %factor_I = 1./apertureInfo.motionModel.indices.nSubPhasePerPosPhase(phase_F);
            
            % store sums of doses (over initial and final phases)
            %dRawCell_sumF{(i-1).*numPhases+phase_I} = dRawCell_sumF{(i-1).*numPhases+phase_I}+dRawTemp.*factor_F;
            %dRawCell_sumI{(i-1).*numPhases+phase_F} = dRawCell_sumI{(i-1).*numPhases+phase_F}+dRawTemp.*factor_I;
            
            for subPhase_F = find(apertureInfo.motionModel.indices.subPhase2PosPhase == phase_F)'
                
                % determine factors
                factor_I = apertureInfo.probI_Ij{i}(phase_I,subPhase_F)./apertureInfo.probI_IJ{i}(phase_I,phase_F);
                
                % store sums of doses (over initial and final phases)
                dRawCell_sumI{(i-1).*numSubPhases+subPhase_F} = dRawCell_sumI{(i-1).*numSubPhases+subPhase_F}+dRawTemp.*factor_I;
            end
            
            for subPhase_I = find(apertureInfo.motionModel.indices.subPhase2PosPhase == phase_I)'
                
                % determine factors
                factor_F = apertureInfo.probF_kL{i}(subPhase_I,phase_F)./apertureInfo.probI_IJ{i}(phase_I,phase_F);
                
                % store sums of doses (over initial and final phases)
                dRawCell_sumF{(i-1).*numSubPhases+subPhase_I} = dRawCell_sumF{(i-1).*numSubPhases+subPhase_I}+dRawTemp.*factor_F;
            end
            
            % do d2 for i1 == i2 now
            % determine probability normalization for d2
            pNorm = pNormMat(phase_I,phase_F);
            
            if abs(pNorm) > eps
                
                % combine into a single multiplication
                d2Raw = dRawTemp.*dRawTemp;
                
                % sum terms
                d2 = d2 + pNorm .* d2Raw;
            end
        end
    end
end

% first loop over all beams
for i1 = 1:numel(apertureInfo.beam)
    
    % second loop over all beams
    for i2 = (i1+1):numel(apertureInfo.beam)
        
        % determine probability normalization matrix
        % if we are looking at different beams, we need to normalize by
        % multiplying by (the probability to transition from F1 to I2)
        % divided (by the probability to arrive at I2)
        
        % calculate the transition and arrival times
        transT12    = sum([apertureInfo.beam((i1+1):(i2-1)).time]);
        
        % calculate the probabilities matrix
        [Pij_transT12,~,~,~] = matRad_transAndTProb(transT12,0,apertureInfo.motionModel);
        
        % first loop over final phases
        for subPhase_F1 = 1:numSubPhases
            
            % get pre-calculated raw doses
            dRaw_1 = dRawCell_sumI{(i1-1).*numSubPhases+subPhase_F1};
            
            % second loop over initial phases
            for subPhase_I2 = 1:numSubPhases
                
                pNorm = Pij_transT12(subPhase_F1,subPhase_I2);
                
                % determine if we're doing d2 for this combo
                if abs(pNorm) > eps
                    
                    % get pre-calculated raw doses
                    dRaw_2      = dRawCell_sumF{(i2-1).*numSubPhases+subPhase_I2};
                    
                    % combine into a single multiplication
                    d2Raw = 2.*dRaw_1.*dRaw_2;
                    
                    % sum terms
                    d2 = d2 + pNorm .* d2Raw;
                end
            end
        end
    end
end

% calculate sum of variances
dVarSmall = d2-d.*d;

dVar = zeros(dij.numOfVoxels,1);
dVar(dij.targetVox) = dVarSmall;

end