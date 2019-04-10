function [dVar,dVarGrad] = matRad_doseVariance(apertureInfo,dij)

%% setup

% get dose from global variable
global matRad_global_d;
d = matRad_global_d;

% allocate gradient of dose
dGrad = zeros(numel(dij.targetVox),size(apertureInfo.apertureVector,1));

% allocate raw doses
dRaw    = cell(apertureInfo.numPhases);
dRaw(:) = {zeros(numel(dij.targetVox),numel(apertureInfo.beam))};

% allocate mean of squared dose
d2 = zeros(dij.numOfVoxels,1);

% allocate gradient of mean of squared dose
d2Grad = zeros(numel(dij.targetVox),size(apertureInfo.apertureVector,1));

% allocate time derivative stuff
pNormGrad   = zeros(1,size(apertureInfo.apertureVector,1));
tIx_Vec     = (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)*apertureInfo.numPhases+(1:apertureInfo.totalNumOfShapes);

%% precalculate doses and gradients
% loop over all beams
for i = 1:numel(apertureInfo.beam)
    
    % loop over phases
    for phase = 1:apertureInfo.numPhases
        
        % pre-calculate raw doses
        % first find relevant bixels
        currBixelIx = apertureInfo.beam(i).bixelIndMap(~isnan(apertureInfo.beam(i).bixelIndMap));
        % now extract bixel weights
        w = apertureInfo.arcI.bixelWeights{phase}(currBixelIx);
        % calculate dose
        dRaw{phase}(:,i) = dij.scaleFactor .* dij.physicalDose{phase}(dij.targetVox,currBixelIx) * w;
        
        % calculate gradient of dose
        % d = d + dij.scaleFactor .*
        % dij.physicalDose{phase_I1}(:,currBixelIx1) * apertureInfo.bixelWeights{phase_I1}(currBixelIx1);
        temp = apertureInfo.bixelJApVec{phase}(:,currBixelIx)';
        dGrad = dGrad + dij.scaleFactor .* dij.physicalDose{phase}(dij.targetVox,currBixelIx) * temp;
        
    end
    
end


% first loop over all beams
for i1 = 1:numel(apertureInfo.beam)
    
    % second loop over all beams
    for i2 = 1:numel(apertureInfo.beam)
        
        %% determine probabilities
        % 1) probability to transfer from the end of the earlier to the
        % beginning of the later
        % 2) probability to arrive at the later
        
        % DOES THIS WORK IF i1 == i2?
        
        % determine which i is first and last
        iFirst = min([i1 i2]);
        iLast = max([i1 i2]);
        
        % calculate the transition and arrival times
        transT12    = sum([apertureInfo.beam((iFirst+1):(iLast-1)).time]);
        T2          = sum([apertureInfo.beam(1:(iLast-1)).time]);
        
        % calculate the probabilities
        [Pij_transT12,Pij_transT12_dot,Pi_T2,Pi_T2_dot] = matRad_transAndTProb(apertureInfo.propVMAT.qij,apertureInfo.propVMAT.initProb,transT12,T2);
        
        % first loop over initial phases
        for phase_I1 = 1:apertureInfo.numPhases
            
            % first loop over final phases
            for phase_F1 = 1:apertureInfo.numPhases
                
                % second loop over initial phases
                for phase_I2 = 1:apertureInfo.numPhases
                    
                    % second loop over final phases
                    for phase_F2 = 1:apertureInfo.numPhases
                        
                        %% determine "weight" tensor quantity, to be multiplied into two dij matrices
                        
                        % first find relevant bixels
                        currBixelIx1 = apertureInfo.beam(i1).bixelIndMap(~isnan(apertureInfo.beam(i1).bixelIndMap));
                        currBixelIx2 = apertureInfo.beam(i2).bixelIndMap(~isnan(apertureInfo.beam(i2).bixelIndMap));
                        
                        % now extract bixel weights corresponding to
                        % initial and final phases of the two beam angles
                        w_I1 = apertureInfo.arcI.bixelWeights{phase_I1}(currBixelIx1);
                        w_F1 = apertureInfo.arcF.bixelWeights{phase_F1}(currBixelIx1);
                        w_I2 = apertureInfo.arcI.bixelWeights{phase_I2}(currBixelIx2);
                        w_F2 = apertureInfo.arcF.bixelWeights{phase_F2}(currBixelIx2);
                        
                        % determine probability normalization
                        pNorm = Pij_transT12(phase_F1,phase_I2)./Pi_T2(phase_I2);
                        
                        % the "weight" tensors are given by of all
                        % combinations of weight products multiplied by the
                        % extra probability to observe both trajectories
                        %W_II = (p_transT12./p_T12).*w_I1*w_I2';
                        %W_IF = (p_transT12./p_T12).*w_I1*w_F2';
                        %W_FI = (p_transT12./p_T12).*w_F1*w_I2';
                        %W_FF = (p_transT12./p_T12).*w_F1*w_F2';
                        
                        %% determine mean of squared dose for this combination
                        
                        % multiply W into dij matrices
                        %d2_II = diag(dij.scaleFactor.^2 * dij.physicalDose{phase_I1}(:,currBixelIx1) * W_II * dij.physicalDose{phase_I2}(:,currBixelIx2)');
                        %d2_IF = diag(dij.scaleFactor.^2 * dij.physicalDose{phase_I1}(:,currBixelIx1) * W_IF * dij.physicalDose{phase_F2}(:,currBixelIx2)');
                        %d2_FI = diag(dij.scaleFactor.^2 * dij.physicalDose{phase_F1}(:,currBixelIx1) * W_FI * dij.physicalDose{phase_I2}(:,currBixelIx2)');
                        %d2_FF = diag(dij.scaleFactor.^2 * dij.physicalDose{phase_F1}(:,currBixelIx1) * W_FF * dij.physicalDose{phase_F2}(:,currBixelIx2)');
                        
                        % pre-calculate raw doses
                        dRaw_I1 = dij.scaleFactor .* dij.physicalDose{phase_I1}(dij.targetVox,currBixelIx1) * w_I1;
                        dRaw_F1 = dij.scaleFactor .* dij.physicalDose{phase_F1}(dij.targetVox,currBixelIx1) * w_F1;
                        dRaw_I2 = dij.scaleFactor .* dij.physicalDose{phase_I2}(dij.targetVox,currBixelIx2) * w_I2;
                        dRaw_F2 = dij.scaleFactor .* dij.physicalDose{phase_F2}(dij.targetVox,currBixelIx2) * w_F2;
                        
                        % use bsxfun instead of repmat? other alternatives?
                        
                        % the "weight" tensors are given by of all
                        % combinations of weight products multiplied by the
                        % extra probability to observe both trajectories
                        % skip tensor determination step and go right to
                        % calculating mean squared doses
                        d2_II = pNorm .* dRaw_I1 .* dRaw_I2;
                        d2_IF = pNorm .* dRaw_I1 .* dRaw_F2;
                        d2_FI = pNorm .* dRaw_F1 .* dRaw_I2;
                        d2_FF = pNorm .* dRaw_F1 .* dRaw_F2;
                        
                        d2_II = dij.scaleFactor.^2 .* pNorm .* (dij.physicalDose{phase_I1}(:,currBixelIx1) * w_I1) .* (dij.physicalDose{phase_I2}(:,currBixelIx2) * w_I2);
                        d2_IF = dij.scaleFactor.^2 .* pNorm .* (dij.physicalDose{phase_I1}(:,currBixelIx1) * w_I1) .* (dij.physicalDose{phase_F2}(:,currBixelIx2) * w_F2);
                        d2_FI = dij.scaleFactor.^2 .* pNorm .* (dij.physicalDose{phase_F1}(:,currBixelIx1) * w_F1) .* (dij.physicalDose{phase_I2}(:,currBixelIx2) * w_I2);
                        d2_FF = dij.scaleFactor.^2 .* pNorm .* (dij.physicalDose{phase_F1}(:,currBixelIx1) * w_F1) .* (dij.physicalDose{phase_F2}(:,currBixelIx2) * w_F2);
                        
                        % sum terms
                        d2 = d2+d2_II+d2_IF+d2_FI+d2_FF;
                        
                        %% now figure out gradients
                        
                        % extract bixel jacobians corresponding to
                        % initial and final phases of the two beam angles
                        j_I1 = apertureInfo.arcI.bixelJApVec{phase_I1}(:,currBixelIx1);
                        j_F1 = apertureInfo.arcF.bixelJApVec{phase_F1}(:,currBixelIx1);
                        j_I2 = apertureInfo.arcI.bixelJApVec{phase_I2}(:,currBixelIx2);
                        j_F2 = apertureInfo.arcF.bixelJApVec{phase_F2}(:,currBixelIx2);
                        
                        % calculate gradient of normalizing probability wrt
                        % times, update gradient vector
                        pNormGrad(tIx_Vec) = sum(apertureInfo.propVMAT.jacobT(:,(iFirst+1):(iLast-1)),2).*Pij_transT12_dot(phase_F1,phase_I2)./Pi_T2(phase_I2) - ...
                            sum(apertureInfo.propVMAT.jacobT(:,1:(iLast-1)),2).*Pi_T2_dot(phase_I2).*Pij_transT12(phase_F1,phase_I2)./(Pi_T2(phase_I2).^2);
                        
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
                        
                        % sum terms
                        d2Grad = d2Grad+d2Grad_II+d2Grad_IF+d2Grad_FI+d2Grad_FF;
                        
                    end
                end
            end
        end
    end
end

% calculate variance of the dose
dVar = d2-d.^2;

% calculate gradient of variance
dVarGrad = d2Grad-2.*d.*dGrad;

end