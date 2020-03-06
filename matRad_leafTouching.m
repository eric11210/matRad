function apertureInfo = matRad_leafTouching(apertureInfo)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to improve instances of leaf touching by moving leaves
% from the centre to sweep with the non-touching leaves.
%
% Currently only works with VMAT, add option to work with IMRT (not as
% crucial)
%
% call
%   apertureInfo = matRad_leafTouching(apertureInfo)
%
% input
%   apertureInfo: matRad aperture weight and shape info struct
%
% output
%   apertureInfo: matRad aperture weight and shape info struct
%
% References
%
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

% for each beam, look for instances where the left and right leaves touch
% replace those positions with the mean leaf position over both banks
% this is to ensure smooth leaf motion

for i = 1:numel(apertureInfo.beam)
    if apertureInfo.propVMAT.beam(i).DAOBeam
        
        if apertureInfo.propVMAT.continuousAperture
            % leafPos defined at ends of arc
            if apertureInfo.propVMAT.beam(i).firstDAO
                % beginning of arc
                
                % extract left and right leaf positions
                leftLeafPos     = apertureInfo.beam(i).shape{1}.leftLeafPos_I_DAO;
                rightLeafPos    = apertureInfo.beam(i).shape{1}.rightLeafPos_I_DAO;
                
                % fix leaf touching
                [leftLeafPos,rightLeafPos] = fixLeafTouching(leftLeafPos,rightLeafPos);

                % replace positions in struct
                apertureInfo.beam(i).shape{1}.leftLeafPos_I_DAO     = leftLeafPos;
                apertureInfo.beam(i).shape{1}.rightLeafPos_I_DAO    = rightLeafPos;
            end
            
            % end of arc
            leftLeafPos = apertureInfo.beam(i).shape{1}.leftLeafPos_F_DAO;
            rightLeafPos = apertureInfo.beam(i).shape{1}.rightLeafPos_F_DAO;
            
            % fix leaf touching
            [leftLeafPos,rightLeafPos] = fixLeafTouching(leftLeafPos,rightLeafPos);
            
            % replace positions in struct
            apertureInfo.beam(i).shape{1}.leftLeafPos_F_DAO     = leftLeafPos;
            apertureInfo.beam(i).shape{1}.rightLeafPos_F_DAO    = rightLeafPos;
        else
            % leafPos defined in middle of arc
            leftLeafPos = apertureInfo.beam(i).shape{1}.leftLeafPos;
            rightLeafPos = apertureInfo.beam(i).shape{1}.rightLeafPos;
            
            % fix leaf touching
            [leftLeafPos,rightLeafPos] = fixLeafTouching(leftLeafPos,rightLeafPos);

            % replace positions in struct
            apertureInfo.beam(i).shape{1}.leftLeafPos   = leftLeafPos;
            apertureInfo.beam(i).shape{1}.rightLeafPos  = rightLeafPos;
        end
    end
end

end

function [leftLeafPos,rightLeafPos] = fixLeafTouching(leftLeafPos,rightLeafPos)
% helper function to fix instances of leaf touching

% find where the leaves touch
touchingInd = leftLeafPos == rightLeafPos;

% take the mean of left and right non-touching leaf positions
meanLeafPos = mean([leftLeafPos(~touchingInd); rightLeafPos(~touchingInd)]);

% substitue touching positions with the mean
leftLeafPos(touchingInd)    = meanLeafPos;
rightLeafPos(touchingInd)   = meanLeafPos;

end

