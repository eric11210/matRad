function stf = matRad_Stf4DPost(stf,masterRayPosBEV,masterRayPosBEV_phase1)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate mask of phase 1 rays for FMO
% generate direct aperture deformation vectors for 4D optimization
%
% call
%   stf = matRad_Stf4DPost(stf,coordsX_vox,coordsY_vox,coordsZ_vox,masterRayPosBEV,numVox)
%
% input
%   stf:        matRad steering information struct
%   pln:        matRad plan meta information struct
%
% output
%   stf:        matRad steering information struct
%
% References
%   -
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

fprintf('matRad: 4D optimization post-processing ... ');

for i = 1:numel(stf)
    stf(i).phase1RayMask = false(size(masterRayPosBEV,1),1);
    
    DAD_phase1 = stf(i).DAD{1};
    
    for phase = 1:numel(stf(i).DAD)
        
        DAD_phaseP = stf(i).DAD{phase};
        stf(i).DAD{phase} = zeros(size(masterRayPosBEV));
        
        for ray = 1:stf(i).numOfRays
            
            if phase == 1
                % this is a mask indicating if a given ray intersects with
                % the target in phase 1
                stf(i).phase1RayMask(ray) = ismember(masterRayPosBEV(ray,:),masterRayPosBEV_phase1,'rows');
            end
            
            % if the DAD is defined for a particular ray, then use it
            % otherwise, set the DAD to be 0 (i.e., it maps to the original
            % ray)
            
            % there are two reasons why a DAD would not be defined for a
            % particular ray: the ray doesn't intersect with the target in
            % 1) that phase or
            % 2) that angle
            
            [DADisDefined,DADind] = ismember(masterRayPosBEV(ray,:),DAD_phase1,'rows');
            if DADisDefined
                stf(i).DAD{phase}(ray,:) = DAD_phaseP(DADind,:);
            else
                stf(i).DAD{phase}(ray,:) = masterRayPosBEV(ray,:);
            end
        end
    end
    
    matRad_progress(i,numel(stf));
end
