function apertureInfo = matRad_doDAD(apertureInfo,stf)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform direct aperture deformation
%
% call
%   apertureInfo = matRad_doDAD(apertureInfo)
%
% input
%   resultGUI:          resultGUI struct to which the output data will be
%                       added, if this field is empty resultGUI struct will
%                       be created
%   stf:                matRad steering information struct
%   dij:                matRad's dij matrix
%   numOfLevels:        number of stratification levels
%   visBool:            toggle on/off visualization (optional)
%
% output
%   resultGUI:          matRad result struct containing the new dose cube
%                       as well as the corresponding weights
%
% References
%   [1] https://doi.org/10.1118/1.2374675
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
shapeInd = 1;
for i = 1:numel(apertureInfo.beam)
    
    % leafPairPos_phase1 is the z-position of each leaf (in phase 1 or P)
    leafPairPos_phase1 = apertureInfo.beam(i).leafPairPos;
    
    % these are the deformation vectors for phase 1
    DADx_phase1 = stf(i).DADx{1};
    DADz_phase1 = stf(i).DADz{1};
    
    shape = apertureInfo.beam(i).shape{1};
    apertureInfo.beam(i).shape = cell(apertureInfo.numPhases,1);
    apertureInfo.beam(i).shape(:) = {shape};
    
    for j = 1:apertureInfo.beam(i).numOfShapes
        
        % these are the x-positions of the various leaves in phase 1
        % we will be querying the DAD at these positions to get the leaf
        % positions at phase P
        leftLeafPos_phase1 = shape(j).leftLeafPos;
        rightLeafPos_phase1 = shape(j).rightLeafPos;
        leftLeafPos_I_phase1 = shape(j).leftLeafPos_I;
        rightLeafPos_I_phase1 = shape(j).rightLeafPos_I;
        leftLeafPos_F_phase1 = shape(j).leftLeafPos_F;
        rightLeafPos_F_phase1 = shape(j).rightLeafPos_F;
        
        for phase = 1:apertureInfo.numPhases
            
            % these are the deformation vectors for phase P
            DADx_phaseP = stf(i).DADx{phase};
            DADz_phaseP = stf(i).DADz{phase};
            
            % use griddata to interpolate the new x- and z-positions of the
            % leaves. We want to know the value of the DAD()_phaseP at the phase 1
            % leaf positions. DADx_phase1 and DADz_phase1 provide a grid.
            leftLeafXPos_phaseP = interp2(DADx_phase1,DADz_phase1,DADx_phaseP,leftLeafPos_phase1,leafPairPos_phase1);
            leftLeafZPos_phaseP = interp2(DADx_phase1,DADz_phase1,DADz_phaseP,leftLeafPos_phase1,leafPairPos_phase1);
            
            rightLeafXPos_phaseP = interp2(DADx_phase1,DADz_phase1,DADx_phaseP,rightLeafPos_phase1,leafPairPos_phase1);
            rightLeafZPos_phaseP = interp2(DADx_phase1,DADz_phase1,DADz_phaseP,rightLeafPos_phase1,leafPairPos_phase1);
            
            leftLeafXPos_I_phaseP = interp2(DADx_phase1,DADz_phase1,DADx_phaseP,leftLeafPos_I_phase1,leafPairPos_phase1);
            leftLeafZPos_I_phaseP = interp2(DADx_phase1,DADz_phase1,DADz_phaseP,leftLeafPos_I_phase1,leafPairPos_phase1);
            
            rightLeafXPos_I_phaseP = interp2(DADx_phase1,DADz_phase1,DADx_phaseP,rightLeafPos_I_phase1,leafPairPos_phase1);
            rightLeafZPos_I_phaseP = interp2(DADx_phase1,DADz_phase1,DADz_phaseP,rightLeafPos_I_phase1,leafPairPos_phase1);
            
            leftLeafXPos_F_phaseP = interp2(DADx_phase1,DADz_phase1,DADx_phaseP,leftLeafPos_F_phase1,leafPairPos_phase1);
            leftLeafZPos_F_phaseP = interp2(DADx_phase1,DADz_phase1,DADz_phaseP,leftLeafPos_F_phase1,leafPairPos_phase1);
            
            rightLeafXPos_F_phaseP = interp2(DADx_phase1,DADz_phase1,DADx_phaseP,rightLeafPos_F_phase1,leafPairPos_phase1);
            rightLeafZPos_F_phaseP = interp2(DADx_phase1,DADz_phase1,DADz_phaseP,rightLeafPos_F_phase1,leafPairPos_phase1);
            
            % the previous interpolation gave us the leaf positions at
            % z-positions not necessarily aligned with the actual MLC. Do
            % another interpolation to correct this. Query the x-positions
            % at the leafPairPos_phase1, with the z-positions as the basis
            leftLeafPos_phaseP = interp1(leftLeafZPos_phaseP,leftLeafXPos_phaseP,leafPairPos_phase1);
            rightLeafPos_phaseP = interp1(rightLeafZPos_phaseP,rightLeafXPos_phaseP,leafPairPos_phase1);
            leftLeafPos_I_phaseP = interp1(leftLeafZPos_I_phaseP,leftLeafXPos_I_phaseP,leafPairPos_phase1);
            rightLeafPos_I_phaseP = interp1(rightLeafZPos_I_phaseP,rightLeafXPos_I_phaseP,leafPairPos_phase1);
            leftLeafPos_F_phaseP = interp1(leftLeafZPos_F_phaseP,leftLeafXPos_F_phaseP,leafPairPos_phase1);
            rightLeafPos_F_phaseP = interp1(rightLeafZPos_F_phaseP,rightLeafXPos_F_phaseP,leafPairPos_phase1);
            
            % fix any instances of NaNs. These occur if the new leaf
            % z-positions (which form the basis of the interpolation)
            % do not cover the entire range of leafPairPos_phase1. In other
            % words, some exterior leaves have been mapped to the interior.
            % This can only happen if the target has moved out of range of
            % those leaves, so the fluence in those rows should be 0.
            % Therefore the left and right leaves should be "touching".
            % Let the touching point be in the middle of the left and right
            % leaves from phase 1.
            ind = isnan(leftLeafPos_phaseP);
            leftLeafPos_phaseP(ind) = (leftLeafPos_phase1(ind)+rightLeafPos_phase1(ind))./2-0.5;
            rightLeafPos_phaseP(ind) = leftLeafPos_phaseP(ind)+1;
            ind = isnan(rightLeafPos_phaseP);
            rightLeafPos_phaseP(ind) = (leftLeafPos_phase1(ind)+rightLeafPos_phase1(ind))./2-0.5;
            leftLeafPos_phaseP(ind) = rightLeafPos_phaseP(ind)+1;
            
            ind = isnan(leftLeafPos_I_phaseP);
            leftLeafPos_I_phaseP(ind) = (leftLeafPos_I_phase1(ind)+rightLeafPos_I_phase1(ind))./2-0.5;
            rightLeafPos_I_phaseP(ind) = leftLeafPos_I_phaseP(ind)+1;
            ind = isnan(rightLeafPos_I_phaseP);
            rightLeafPos_I_phaseP(ind) = (leftLeafPos_I_phase1(ind)+rightLeafPos_I_phase1(ind))./2-0.5;
            leftLeafPos_I_phaseP(ind) = rightLeafPos_I_phaseP(ind)+1;
            
            ind = isnan(leftLeafPos_F_phaseP);
            leftLeafPos_F_phaseP(ind) = (leftLeafPos_F_phase1(ind)+rightLeafPos_F_phase1(ind))./2-0.5;
            rightLeafPos_F_phaseP(ind) = leftLeafPos_F_phaseP(ind)+1;
            ind = isnan(rightLeafPos_F_phaseP);
            rightLeafPos_F_phaseP(ind) = (leftLeafPos_F_phase1(ind)+rightLeafPos_F_phase1(ind))./2-0.5;
            leftLeafPos_F_phaseP(ind) = rightLeafPos_F_phaseP(ind)+1;
            
            % now put the new leaf positions in apertureInfo
            apertureInfo.beam(i).shape{phase}(j).leftLeafPos = leftLeafPos_phaseP;
            apertureInfo.beam(i).shape{phase}(j).rightLeafPos = rightLeafPos_phaseP;
            apertureInfo.beam(i).shape{phase}(j).leftLeafPos_I = leftLeafPos_I_phaseP;
            apertureInfo.beam(i).shape{phase}(j).rightLeafPos_I = rightLeafPos_I_phaseP;
            apertureInfo.beam(i).shape{phase}(j).leftLeafPos_F = leftLeafPos_F_phaseP;
            apertureInfo.beam(i).shape{phase}(j).rightLeafPos_F = rightLeafPos_F_phaseP;
            
            % fix the vectorOffsets
            % there are two shifts: one to make room for the weights, and
            % another to make room for the new leaf positions at each
            % phase
            apertureInfo.beam(i).shape{phase}(j).vectorOffset = apertureInfo.beam(i).shape{phase}(j).vectorOffset + (apertureInfo.numPhases-1)*apertureInfo.totalNumOfShapes + (phase-1)*apertureInfo.totalNumOfLeafPairs;
            apertureInfo.beam(i).shape{phase}(j).weightOffset = apertureInfo.beam(i).shape{phase}(j).weightOffset + (phase-1)*apertureInfo.totalNumOfShapes;
        end
        
        %% TEMP
        apertureInfo.propVMAT.beam(i).timeInd = apertureInfo.numPhases*(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)+shapeInd;
        shapeInd = shapeInd+1;
    end
end

% now, initial and final leaf positions of neighbouring angles may not be
% equal (different DAD at different angles)
% loop over phases and angles to fix these instances (reset both I and F to
% be the mean of I and F).
for i = 1:numel(apertureInfo.beam)
    for phase = 1:apertureInfo.numPhases
        if apertureInfo.propVMAT.continuousAperture
            if apertureInfo.propVMAT.beam(i).DAOBeam && ~apertureInfo.propVMAT.beam(i).doseAngleDAO(2)
                % if this is true, then the F of beam i should be equal to the I of
                % beam i+1
                % set them both to the mean
                temp_leftLeafPos_F = apertureInfo.beam(i).shape{phase}(j).leftLeafPos_F;
                temp_rightLeafPos_F = apertureInfo.beam(i).shape{phase}(j).rightLeafPos_F;
                temp_leftLeafPos_I = apertureInfo.beam(i+1).shape{phase}(j).leftLeafPos_I;
                temp_rightLeafPos_I = apertureInfo.beam(i+1).shape{phase}(j).rightLeafPos_I;
                
                new_leftLeafPos = mean([temp_leftLeafPos_I, temp_leftLeafPos_F],2);
                new_rightLeafPos = mean([temp_rightLeafPos_I, temp_rightLeafPos_F],2);
                
                apertureInfo.beam(i).shape{phase}(j).leftLeafPos_F = new_leftLeafPos;
                apertureInfo.beam(i).shape{phase}(j).rightLeafPos_F = new_rightLeafPos;
                apertureInfo.beam(i+1).shape{phase}(j).leftLeafPos_I = new_leftLeafPos;
                apertureInfo.beam(i+1).shape{phase}(j).rightLeafPos_I = new_rightLeafPos;
                
            end
        end
        
        for j = 1:apertureInfo.beam(i).numOfShapes
            % fix any instances where the right leaf is somehow to the left
            % of the left leaf. Set them both to be the mean
            
            leftLeafPos = apertureInfo.beam(i).shape{phase}(j).leftLeafPos;
            rightLeafPos = apertureInfo.beam(i).shape{phase}(j).rightLeafPos;
            leftLeafPos_I = apertureInfo.beam(i).shape{phase}(j).leftLeafPos_I;
            rightLeafPos_I = apertureInfo.beam(i).shape{phase}(j).rightLeafPos_I;
            leftLeafPos_F = apertureInfo.beam(i).shape{phase}(j).leftLeafPos_F;
            rightLeafPos_F = apertureInfo.beam(i).shape{phase}(j).rightLeafPos_F;
            
            ind = rightLeafPos < leftLeafPos;
            apertureInfo.beam(i).shape{phase}(j).leftLeafPos(ind) = mean([leftLeafPos(ind), rightLeafPos(ind)],2);
            apertureInfo.beam(i).shape{phase}(j).rightLeafPos(ind) = mean([leftLeafPos(ind), rightLeafPos(ind)],2);
            
            ind = rightLeafPos_I < leftLeafPos_I;
            apertureInfo.beam(i).shape{phase}(j).leftLeafPos_I(ind) = mean([leftLeafPos_I(ind), rightLeafPos_I(ind)],2);
            apertureInfo.beam(i).shape{phase}(j).rightLeafPos_I(ind) = mean([leftLeafPos_I(ind), rightLeafPos_I(ind)],2);
            
            ind = rightLeafPos_F < leftLeafPos_F;
            apertureInfo.beam(i).shape{phase}(j).leftLeafPos_F(ind) = mean([leftLeafPos_F(ind), rightLeafPos_F(ind)],2);
            apertureInfo.beam(i).shape{phase}(j).rightLeafPos_F(ind) = mean([leftLeafPos_F(ind), rightLeafPos_F(ind)],2);
            
        end
    end
end
