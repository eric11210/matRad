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

fprintf('matRad: 4D optimization post-processing ... ');

for i = 1:numel(stf)
    stf(i).phase1RayMask = false(size(masterRayPosBEV,1),1);
    stf(i).DADx = cell(numel(stf(i).DAD),1);
    stf(i).DADz = cell(numel(stf(i).DAD),1);
    
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
        
        % extend the domain of the DAD so that it is a rectangle covering
        % all rays, plus one more
        % also change the format of the DAD to be a grid
        minX = min(stf(1).DAD{1}(:,1))-stf(i).bixelWidth;
        maxX = max(stf(1).DAD{1}(:,1))+stf(i).bixelWidth;
        minZ = min(stf(1).DAD{1}(:,3))-stf(i).bixelWidth;
        maxZ = max(stf(1).DAD{1}(:,3))+stf(i).bixelWidth;
        X = minX:stf(i).bixelWidth:maxX;
        Z = minZ:stf(i).bixelWidth:maxZ;
        
        sizX = (maxX-minX)./stf(i).bixelWidth+1;
        sizZ = (maxZ-minZ)./stf(i).bixelWidth+1;
        
        stf(i).DADx{phase} = zeros(sizZ,sizX);
        stf(i).DADz{phase} = zeros(sizZ,sizX);
        
        for zz = 1:sizZ
            for xx = 1:sizX
                
                x = X(xx);
                z = Z(zz);
                [DADisDefined,DADind] = ismember([x 0 z],stf(i).DAD{1},'rows');
                if DADisDefined
                    stf(i).DADx{phase}(zz,xx) = stf(i).DAD{phase}(DADind,1);
                    stf(i).DADz{phase}(zz,xx) = stf(i).DAD{phase}(DADind,3);
                else
                    stf(i).DADx{phase}(zz,xx) = x;
                    stf(i).DADz{phase}(zz,xx) = z;
                end
                
            end
        end
    end
    
    matRad_progress(i,numel(stf));
end
