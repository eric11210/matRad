function apertureInfo = matRad_apertures2Library(apertureInfo,pln,stf,numPhases)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert a sequence of apertures to a full library, using the same
% aperture for each phase
%
% call
%   apertureInfo = matRad_sequence2Library(apertureInfo,numPhases)
%
% input
%   apertureInfo:       matRad aperture information struct
%   pln:                matRad plan struct
%   numPhases:          number of phases for the library
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

apertureInfo.numPhases  = numPhases;
pln.propOpt.run4D       = true;

% determine if all phases of the library already exist
% this means that we have previously converted the apertures to a library
allPhases = numel(apertureInfo.beam(1).shape) == apertureInfo.numPhases;

if ~allPhases
    
    shapeInd = 1;
    
    for i = 1:numel(apertureInfo.beam)
        
        % copy shape properties, including leaf positions and weight
        shape = apertureInfo.beam(i).shape{1};
        apertureInfo.beam(i).shape = cell(apertureInfo.numPhases,1);
        apertureInfo.beam(i).shape(:) = {shape};
        
        if apertureInfo.propVMAT.beam(i).DAOBeam
            % fix the time indices
            apertureInfo.propVMAT.beam(i).timeInd = apertureInfo.numPhases*(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)+shapeInd;
            shapeInd = shapeInd+1;
            
            for j = 1:apertureInfo.beam(i).numOfShapes
                for phase = 1:apertureInfo.numPhases
                    
                    % fix the vectorOffsets
                    % there are two shifts: one to make room for the weights, and
                    % another to make room for the new leaf positions at each
                    % phase
                    apertureInfo.beam(i).shape{phase}(j).vectorOffset = apertureInfo.beam(i).shape{phase}(j).vectorOffset + (apertureInfo.numPhases-1)*apertureInfo.totalNumOfShapes + (phase-1)*apertureInfo.totalNumOfLeafPairs;
                    
                    % also fix the weightOffsets
                    % there is one shift, to make room for the weights
                    apertureInfo.beam(i).shape{phase}(j).weightOffset = apertureInfo.beam(i).shape{phase}(j).weightOffset + (phase-1)*apertureInfo.totalNumOfShapes;
                end
            end
        end
    end
    
    % refresh metadata
    apertureInfo = matRad_apertureInfoMeta(apertureInfo,pln,stf,false);
    
    % update vector
    [apertureInfo.apertureVector, apertureInfo.mappingMx, apertureInfo.limMx] = matRad_daoApertureInfo2Vec(apertureInfo);
end

% update  apertureInfo
apertureInfo = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfo.apertureVector);

end