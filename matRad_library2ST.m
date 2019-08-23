function apertureInfo = matRad_library2ST(apertureInfo,pln,stf,trajectory)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract a sequence of apertures from a full library using the breathing
% trajectory specified in the trajectory structure
%
% call
%   apertureInfo = matRad_library2ST(apertureInfo,trajectory)
%
% input
%   apertureInfo:       matRad aperture information struct
%   trajectory:         trajectory information struct
%
% output
%   apertureInfo:       matRad aperture information struct
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

numPhasesReal = apertureInfo.numPhases;

apertureInfo.run4D          = false;
apertureInfo.numPhases      = 1;
apertureInfo.jacobiScale    = ones(apertureInfo.totalNumOfShapes.*apertureInfo.numPhases,1);

shapeInd = 1;

for i = 1:numel(apertureInfo.beam)
    
    % determine correct phase for current angle
    phase = trajectory.pSimulated(i);
    
    % copy shape properties, including leaf positions and weight
    shape = apertureInfo.beam(i).shape{phase};
    apertureInfo.beam(i).shape = cell(1,1);
    apertureInfo.beam(i).shape{1} = shape;
    
    if apertureInfo.propVMAT.beam(i).DAOBeam
        % fix the time indices
        apertureInfo.propVMAT.beam(i).timeInd = apertureInfo.numPhases*(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)+shapeInd;
        shapeInd = shapeInd+1;
        
        for j = 1:apertureInfo.beam(i).numOfShapes
            
            % fix the vectorOffsets
            % there are two shifts: one to compensating for the lost weights,
            % and another to comensate for the lost leaf positions at each
            % phase
            apertureInfo.beam(i).shape{1}(j).vectorOffset = shape(j).vectorOffset - (numPhasesReal-1)*apertureInfo.totalNumOfShapes - (phase-1)*apertureInfo.totalNumOfLeafPairs;
            
            % also fix the weightOffsets
            % there is one shift, to compensate for the lost weights
            apertureInfo.beam(i).shape{1}(j).weightOffset = shape(j).weightOffset - (phase-1)*apertureInfo.totalNumOfShapes;
            
            % correct jacobi scale
            apertureInfo.jacobiScale(apertureInfo.beam(i).shape{1}(j).weightOffset) = apertureInfo.beam(i).shape{1}(j).jacobiScale;
        end
    end
end

% refresh metadata
apertureInfo = matRad_apertureInfoMeta(apertureInfo,pln,stf,apertureInfo.totalNumOfOptBixels);

% update vector
[apertureInfo.apertureVector, apertureInfo.mappingMx, apertureInfo.limMx] = matRad_daoApertureInfo2Vec(apertureInfo);

% update  apertureInfo
apertureInfo = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfo.apertureVector);

end