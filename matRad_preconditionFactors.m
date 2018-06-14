function apertureInfo = matRad_preconditionFactors(apertureInfo)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate preconditioning factors for DAO (only what Esther Wild called
% the Jacobi preconditioner, not the dij scaled). Scale weights in vector
% accordingly.
%
% call
%   apertureInfo =
%   matRad_preconditionFactors(apertureInfo)
%
% input
%   apertureInfo:       aperture shape info struct
%
% output
%   apertureInfo:       aperture shape info struct with new factors
%
% References
%   [1] http://onlinelibrary.wiley.com/doi/10.1118/1.4914863/full
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

% This is the dij scaling factor which will be applied during DAO.  It is
% given by the dividing the mean of the actual aperture weights by the
% bixel width.  This factor will divide all of the aperture weights.
if ~apertureInfo.propVMAT.continuousAperture
    dijScaleFactor = mean(apertureInfo.apertureVector(1:apertureInfo.totalNumOfShapes*apertureInfo.numPhases))/(apertureInfo.bixelWidth);
else
    dijScaleFactor = mean(apertureInfo.apertureVector(1:apertureInfo.totalNumOfShapes*apertureInfo.numPhases))/(2*apertureInfo.bixelWidth);
end

apertureInfo.jacobiScale = zeros(apertureInfo.totalNumOfShapes*apertureInfo.numPhases,1);

shapeInd = 1;

for phase = 1:apertureInfo.numPhases
    for i = 1:numel(apertureInfo.beam)
        
        if ~apertureInfo.runVMAT || (apertureInfo.runVMAT && apertureInfo.propVMAT.beam(i).DAOBeam)
            % in other words, do this for every beam if it's not VMAT, and for
            % optimized beams only if it is
            
            for j = 1:apertureInfo.beam(i).numOfShapes
                
                % To get the jacobi scaling factor, first factor the
                % current aperture's weight out of the dijScaling factor.  Also
                % remove the bixel width.  Now we have the mean weight relative
                % to the current weight.
                % Next, multiply by the sqrt of ~approximately the number of
                % open bixels (slight modification to Esther Wild's formula).
                % The variables corresponding to the aperture weights will be
                % multiplied by this number, which will decrease the gradients.
                if ~apertureInfo.propVMAT.continuousAperture
                    apertureInfo.beam(i).shape{phase}(j).jacobiScale = (dijScaleFactor.*apertureInfo.bixelWidth./apertureInfo.beam(i).shape{phase}(j).weight).*sqrt(sum(apertureInfo.beam(i).shape{phase}(j).shapeMap(:).^2));
                else
                    apertureInfo.beam(i).shape{phase}(j).jacobiScale = (dijScaleFactor./apertureInfo.beam(i).shape{phase}(j).weight).*sqrt(sum(apertureInfo.beam(i).shape{phase}(j).shapeMap(:).^2))./apertureInfo.beam(i).shape{phase}(j).sqrtSumGradSq;
                end
                apertureInfo.jacobiScale(shapeInd) = apertureInfo.beam(i).shape{phase}(j).jacobiScale;
                
                apertureInfo.apertureVector(shapeInd) = apertureInfo.beam(i).shape{phase}(j).jacobiScale*apertureInfo.beam(i).shape{phase}(j).weight;
                
                shapeInd = shapeInd+1;
            end
        end
    end
end



