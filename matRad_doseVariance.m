function dVar = matRad_doseVariance(apertureInfo,dij)

%% setup

% allocate raw doses
dRaw = cell(numel(apertureInfo.beam).*apertureInfo.numPhases.^2,1);

% allocate mean of dose
d = zeros(numel(dij.structVox),1);

% allocate mean of squared dose
d2 = zeros(numel(dij.structVox),1);

%% precalculate doses and gradients
% loop over all beams
for i = 1:numel(apertureInfo.beam)
    
    % pre-calculate raw doses
    % first find relevant bixels
    currBixelIx = apertureInfo.beam(i).bixelIndMap(~isnan(apertureInfo.beam(i).bixelIndMap));
    
    % loop over initial and final phases
    for phase_I = 1:apertureInfo.numPhases
        
        % extract dijBeamPhase for phase_I
        dijBeamPhase_I      = dij.scaleFactor .* dij.physicalDose{phase_I}(dij.structVox,currBixelIx);
        
        for phase_F = 1:apertureInfo.numPhases
            
            % determine cell index
            cellInd = (phase_I-1).*apertureInfo.numPhases+phase_F;
            
            % extract dijBeamPhase for phase_F
            if phase_I == phase_F
                dijBeamPhase_F      = dijBeamPhase_I;
            else
                dijBeamPhase_F      = dij.scaleFactor .* dij.physicalDose{phase_F}(dij.structVox,currBixelIx);
            end
            
            % now extract bixel weights
            w_I = apertureInfo.arcI.bixelWeights{cellInd}(currBixelIx);
            w_F = apertureInfo.arcF.bixelWeights{cellInd}(currBixelIx);
            % calculate dose
            dRawTemp_I  = dijBeamPhase_I * w_I;
            dRawTemp_F  = dijBeamPhase_F * w_F;
            dRawTemp    = dRawTemp_I+dRawTemp_F;
            
            % accumulate sum of dose
            d = d+dRawTemp;
            
            % store doses
            dRaw{(i-1).*apertureInfo.numPhases.^2+cellInd} = dRawTemp;
            
        end
    end
end

% first loop over all beams
for i1_loop = 1:numel(apertureInfo.beam)
    
    % second loop over all beams
    for i2_loop = 1:numel(apertureInfo.beam)
        
        %% housekeeping
        
        % determine the real i1 and i2 from the loop variables
        % hint: i1 must be smaller than or equal to i2
        i1 = min([i1_loop i2_loop]);
        i2 = max([i1_loop i2_loop]);
        
        %% determine probabilities
        
        % determine probability normalization matrix
        if i1 ~= i2
            % if we are looking at different beams, we need to normalize by
            % multiplying by (the probability to transition from F1 to I2)
            % divided (by the probability to arrive at I2)
            
            % calculate the transition and arrival times
            transT12    = sum([apertureInfo.beam((i1+1):(i2-1)).time]);
            T2          = sum([apertureInfo.beam(1:(i2-1)).time]);
            
            % calculate the probabilities
            [Pij_transT12,~,Pi_T2,~] = matRad_transAndTProb(transT12,T2,apertureInfo.motionModel);
            
            % calculate the matrix
            %pNorm      = Pij_transT12(phase_F1,phase_I2)./Pi_T2(phase_I2);
            pNormMat    = Pij_transT12./Pi_T2;
        else
            % if we are looking at the same beam, we need to normalize by
            % dividing by (the probability to arrive at I1) multiplied by
            % (the probability to transition from I1 to F1)
            
            % calculate the transition and arrival times
            transT  = apertureInfo.beam(i).time;
            T       = sum([apertureInfo.beam(1:(i-1)).time]);
            
            % calculate the transition and arrival times
            [Pij_transT,~,Pi_T,~] = matRad_transAndTProb(transT,T,apertureInfo.motionModel);
            
            % calculate the matrix
            %pNorm      = 1./(Pi_T(phase_I1).*Pij_transT(phase_I1,phase_F1));
            pNormMat    = 1./(Pi_T'.*Pij_transT);
        end
        
        % first loop over initial phases
        for phase_I1 = 1:apertureInfo.numPhases
            
            % first loop over final phases
            for phase_F1 = 1:apertureInfo.numPhases
                
                % determine cell index for beam 1
                cellInd_1 = (phase_I1-1).*apertureInfo.numPhases+phase_F1;
                
                % get pre-calculated raw doses
                dRaw_1      = dRaw{(i1-1).*apertureInfo.numPhases.^2+cellInd_1};
                
                if i1 ~= i2
                    % if we are looking at different beams, then we need to
                    % loop over initial and final phases for the second
                    % beam
                    
                    % second loop over initial phases
                    for phase_I2 = 1:apertureInfo.numPhases
                        
                        % determine probability normalization
                        pNorm = pNormMat(phase_F1,phase_I2);
                        
                        % determine if we're doing d2 for this combo
                        dod2 = abs(pNorm) > eps;
                        
                        % second loop over final phases
                        for phase_F2 = 1:apertureInfo.numPhases
                            
                            % determine cell index for beam 2
                            cellInd_2 = (phase_I2-1).*apertureInfo.numPhases+phase_F2;
                            
                            %% determine mean of squared dose for this combination
                            % only do the following if we need d2 stuff,
                            % i.e. if pNorm is > eps
                            if dod2
                                
                                % get pre-calculated raw doses
                                dRaw_2      = dRaw{(i2-1).*apertureInfo.numPhases.^2+cellInd_2};
                                
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
                                d2SumRaw = dRaw_1.*dRaw_2;
                                
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
                                d2 = d2 + pNorm .* d2SumRaw;
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
                else
                    % if we are looking at the same beam, i1 == i2, then
                    % the initial and final phases for the second beam are
                    % necessarily the same as the first
                    
                    % determine probability normalization
                    pNorm = pNormMat(phase_I1,phase_F1);
                    
                    % determine if we're doing d2 for this combo
                    dod2 = abs(pNorm) > eps;
                    
                    %% determine mean of squared dose for this combination
                    % only do the following if we need d2 stuff,
                    % i.e. if pNorm is > eps
                    if dod2
                        % here, dRaw2 = dRaw1
                        
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
                        d2SumRaw = dRaw_1.*dRaw_1;
                        
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
                        d2 = d2 + pNorm .* d2SumRaw;
                    end
                end
            end
        end
    end
end

% calculate sum of variances
dVar = d2-d.*d;

end