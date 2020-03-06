function apertureInfo = matRad_optDelivery(apertureInfo,fast)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad: optimize VMAT delivery
%
% call
%   matRad_optDelivery(result,fast)
%
% input
%   result:             result struct from fluence optimization/sequencing
%   fast:               1 => fastest possible delivery
%                       0 => mutliply delivery time by 10%
%
% output
%   apertureInfo:       aperture shape info struct
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    fast = 1;
end


%speed up delivery time, when it is permitted by constraints
%constraints to consider: doserate, leaf speed, and gantry speed

%Do this after DAO

fileName = apertureInfo.propVMAT.machineConstraintFile;
try
    load(fileName,'machine');
catch
    error(['Could not find the following machine file: ' fileName ]);
end

%calculate max leaf speed
apertureInfo = matRad_maxLeafSpeed(apertureInfo);

for i = 1:size(apertureInfo.beam,2)
    
    if ~apertureInfo.propVMAT.fixedGantrySpeed
        % only change times if we're not using a fixed gantry speed
        
        %all of these should be greater than 1, since DAO respects the
        %constraints
        
        %if one of them is less than 1, then a constraint is violated
        factorMURate = inf;
        for phase = 1:apertureInfo.numPhases
            factorMURate = min([machine.constraints.monitorUnitRate(2)/apertureInfo.beam(i).shape{phase}.MURate factorMURate]);
        end
        factorLeafSpeed = machine.constraints.leafSpeed(2)/apertureInfo.beam(i).maxLeafSpeed;
        factorGantryRot = machine.constraints.gantryRotationSpeed(2)/apertureInfo.beam(i).gantryRot;
        
        %The constraint that is limiting the speed the most is the one
        %whose factor is closest to 1
        factor = min([factorMURate factorLeafSpeed factorGantryRot]);
        if ~fast
            %if the limiting rate is already 10% lower than the limit,
            %then do nothing (factor = 1)
            %otherwise, scale rates so that the limiting rate is 10% lower
            %than the limit
            factor = min([1 factor*0.9]);
        end
        
        %multiply each speed by this factor
        for phase = 1:apertureInfo.numPhases
            apertureInfo.beam(i).shape{phase}.MURate = factor*apertureInfo.beam(i).shape{phase}.MURate;
        end
        apertureInfo.beam(i).maxLeafSpeed = factor*apertureInfo.beam(i).maxLeafSpeed;
        apertureInfo.beam(i).gantryRot = factor*apertureInfo.beam(i).gantryRot;
        apertureInfo.beam(i).time = apertureInfo.beam(i).time/factor;
    end
    
    for phase = 1:apertureInfo.numPhases
        factorsMURate = machine.constraints.monitorUnitRate/apertureInfo.beam(i).shape{phase}.MURate;
        
        if factorsMURate(1) > 1
            apertureInfo.beam(i).shape{phase}.MURate = factorsMURate(1)*apertureInfo.beam(i).shape{phase}.MURate;
            apertureInfo.beam(i).shape{phase}(1).weight = factorsMURate(1)*apertureInfo.beam(i).shape{phase}(1).weight;
        end
        
        if factorsMURate(2) < 1
            apertureInfo.beam(i).shape{phase}.MURate = factorsMURate(2)*apertureInfo.beam(i).shape{phase}.MURate;
            apertureInfo.beam(i).shape{phase}(1).weight = factorsMURate(2)*apertureInfo.beam(i).shape{phase}(1).weight;
        end
    end
end

%recalculate vector with new times and weights
[apertureInfo.apertureVector,~,~] = matRad_daoApertureInfo2Vec(apertureInfo);

%redo bixel weight calculation
apertureInfo = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfo.apertureVector);

end

