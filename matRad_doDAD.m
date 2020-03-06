function apertureInfo = matRad_doDAD(apertureInfo,stf)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform direct aperture deformation
%
% call
%   apertureInfo = matRad_doDAD(apertureInfo,stf)
%
% input
%   apertureInfo:       matRad aperture information struct
%   stf:                matRad steering information struct
%
% output
%   apertureInfo:       matRad aperture information struct
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

% determine if all phases of the library already exist
% this means that we have previously done a DAD-like transformation
allPhases = numel(apertureInfo.beam(1).shape) == apertureInfo.numPhases;

for i = 1:numel(apertureInfo.beam)
    
    if apertureInfo.propVMAT.beam(i).DAOBeam
        
        % leafPairPos_phase1 is the z-position of each leaf (in phase 1 or P)
        leafPairPos_phase1 = apertureInfo.beam(i).leafPairPos;
        
        % these are the deformation vectors for phase 1
        DADx_phase1 = stf(i).DADx{1};
        DADz_phase1 = stf(i).DADz{1};
        
        shape = apertureInfo.beam(i).shape{1};
        
        if ~allPhases && numel(apertureInfo.beam(i).shape) == 1
            % if all the phases don't exist already, create them
            apertureInfo.beam(i).shape = cell(apertureInfo.numPhases,1);
            apertureInfo.beam(i).shape(:) = {shape};
        elseif ~allPhases
            % throw error
            error('Number of phases in library is more than 1, but not all phases exist');
        end
        
        for j = 1:numel(shape)
            
            for phase = 1:apertureInfo.numPhases
                
                % these are the deformation vectors for phase P
                DADx_phaseP = stf(i).DADx{phase};
                DADz_phaseP = stf(i).DADz{phase};
                
                if apertureInfo.propVMAT.continuousAperture
                    % leafPos defined at ends of arc
                    if apertureInfo.propVMAT.beam(i).firstDAO
                        % beginning of arc
                        
                        % extract left and right leaf positions
                        leftLeafPos_phase1  = shape(j).leftLeafPos_I_DAO;
                        rightLeafPos_phase1 = shape(j).rightLeafPos_I_DAO;
                        
                        % do DAD
                        [leftLeafPos_phaseP,rightLeafPos_phaseP] = doDAD(leftLeafPos_phase1,rightLeafPos_phase1,leafPairPos_phase1,DADx_phase1,DADz_phase1,DADx_phaseP,DADz_phaseP);
                        
                        % replace positions in struct
                        apertureInfo.beam(i).shape{phase}(j).leftLeafPos_I_DAO  = leftLeafPos_phaseP;
                        apertureInfo.beam(i).shape{phase}(j).rightLeafPos_I_DAO = rightLeafPos_phaseP;
                    end
                    
                    % end of arc
                    leftLeafPos_phase1  = shape(j).leftLeafPos_F_DAO;
                    rightLeafPos_phase1 = shape(j).rightLeafPos_F_DAO;
                    
                    % do DAD
                    [leftLeafPos_phaseP,rightLeafPos_phaseP] = doDAD(leftLeafPos_phase1,rightLeafPos_phase1,leafPairPos_phase1,DADx_phase1,DADz_phase1,DADx_phaseP,DADz_phaseP);
                    
                    % replace positions in struct
                    apertureInfo.beam(i).shape{phase}(j).leftLeafPos_F_DAO  = leftLeafPos_phaseP;
                    apertureInfo.beam(i).shape{phase}(j).rightLeafPos_F_DAO = rightLeafPos_phaseP;
                else
                    % leafPos defined in middle of arc
                    leftLeafPos_phase1  = shape(j).leftLeafPos;
                    rightLeafPos_phase1 = shape(j).rightLeafPos;
                    
                    % do DAD
                    [leftLeafPos_phaseP,rightLeafPos_phaseP] = doDAD(leftLeafPos_phase1,rightLeafPos_phase1,leafPairPos_phase1,DADx_phase1,DADz_phase1,DADx_phaseP,DADz_phaseP);
                    
                    % replace positions in struct
                    apertureInfo.beam(i).shape{phase}(j).leftLeafPos    = leftLeafPos_phaseP;
                    apertureInfo.beam(i).shape{phase}(j).rightLeafPos   = rightLeafPos_phaseP;
                end
                
                if ~allPhases && apertureInfo.propVMAT.beam(i).DAOBeam
                    % fix the vectorOffsets
                    % there are two shifts: one to make room for the weights, and
                    % another to make room for the new leaf positions at each
                    % phase
                    apertureInfo.beam(i).shape{phase}(j).vectorOffset = apertureInfo.beam(i).shape{phase}(j).vectorOffset + (apertureInfo.numPhases-1)*apertureInfo.totalNumOfShapes + (phase-1)*apertureInfo.totalNumOfLeafPairs;
                    apertureInfo.beam(i).shape{phase}(j).weightOffset = apertureInfo.beam(i).shape{phase}(j).weightOffset + (phase-1)*apertureInfo.totalNumOfShapes;
                end
            end
        end
    end
end

end


function [leftLeafPos_phaseP,rightLeafPos_phaseP] = doDAD(leftLeafPos_phase1,rightLeafPos_phase1,leafPairPos_phase1,DADx_phase1,DADz_phase1,DADx_phaseP,DADz_phaseP)

% use griddata to interpolate the new x- and z-positions of the
% leaves. We want to know the value of the DAD()_phaseP at the phase 1
% leaf positions. DADx_phase1 and DADz_phase1 provide a grid.
leftLeafXPos_phaseP     = griddata(DADx_phase1,DADz_phase1,DADx_phaseP,leftLeafPos_phase1,leafPairPos_phase1);
leftLeafZPos_phaseP     = griddata(DADx_phase1,DADz_phase1,DADz_phaseP,leftLeafPos_phase1,leafPairPos_phase1);
rightLeafXPos_phaseP    = griddata(DADx_phase1,DADz_phase1,DADx_phaseP,rightLeafPos_phase1,leafPairPos_phase1);
rightLeafZPos_phaseP    = griddata(DADx_phase1,DADz_phase1,DADz_phaseP,rightLeafPos_phase1,leafPairPos_phase1);

% the previous interpolation gave us the leaf positions at
% z-positions not necessarily aligned with the actual MLC. Do
% another interpolation to correct this. Query the x-positions
% at the leafPairPos_phase1, with the z-positions as the basis
leftLeafPos_phaseP  = interp1(leftLeafZPos_phaseP,leftLeafXPos_phaseP,leafPairPos_phase1);
rightLeafPos_phaseP = interp1(rightLeafZPos_phaseP,rightLeafXPos_phaseP,leafPairPos_phase1);

% fix any instances of NaNs. These occur if the new leaf
% z-positions (which form the basis of the interpolation)
% do not cover the entire range of leafPairPos_phase1. In other
% words, some exterior leaves have been mapped to the interior.
% This can only happen if the target has moved out of range of
% those leaves, so the fluence in those rows should be 0.
% Therefore the left and right leaves should be "touching".
% also identify any leaves which are actually still touching
touchingInd = isnan(leftLeafPos_phaseP) | isnan(rightLeafPos_phaseP)| leftLeafPos_phaseP == rightLeafPos_phaseP | rightLeafPos_phaseP-leftLeafPos_phaseP < 1e-8;

% Let the touching point be the mean of the non-touching leaf positions
meanLeafPos = mean([leftLeafPos_phaseP(~touchingInd); rightLeafPos_phaseP(~touchingInd)]);

% substitue touching positions with the mean
leftLeafPos_phaseP(touchingInd)    = meanLeafPos;
rightLeafPos_phaseP(touchingInd)   = meanLeafPos;

end