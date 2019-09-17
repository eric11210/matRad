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
    
    for phase = 1:numel(stf(i).DAD)
        
        % redefine the DADs to be in the same order as the masterRayPosBEV
        DAD_phaseP          = stf(i).DAD{phase};
        stf(i).DAD{phase}   = zeros(size(masterRayPosBEV));
        
        % loop through all the rays
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
            
            [DADisDefined,DADind] = ismember(masterRayPosBEV(ray,:),stf(i).rayPos_phase1,'rows');
            if DADisDefined
                stf(i).DAD{phase}(ray,:) = DAD_phaseP(DADind,:);
            else
                stf(i).DAD{phase}(ray,:) = NaN;
            end
        end
        
        % extend the domain of the DAD so that it is a rectangle covering
        % all rays, plus one more
        % also change the format of the DAD to be a grid
        minX = min(masterRayPosBEV(:,1))-stf(i).bixelWidth;
        maxX = max(masterRayPosBEV(:,1))+stf(i).bixelWidth;
        minZ = min(masterRayPosBEV(:,3))-stf(i).bixelWidth;
        maxZ = max(masterRayPosBEV(:,3))+stf(i).bixelWidth;
        X = minX:stf(i).bixelWidth:maxX;
        Z = minZ:stf(i).bixelWidth:maxZ;
        
        sizX = (maxX-minX)./stf(i).bixelWidth+1;
        sizZ = (maxZ-minZ)./stf(i).bixelWidth+1;
        
        % we want the DADx/z matrices to transition smoothly from their
        % correct values in the target to identity at the edges
        
        % define the DAD arrays
        stf(i).DADx{phase} = zeros(sizZ,sizX);
        stf(i).DADz{phase} = zeros(sizZ,sizX);
        
        % define the x and z positions of the DADx/z_phase0
        DADx_phase0 = zeros(sizZ,sizX);
        DADz_phase0 = zeros(sizZ,sizX);
        
        % this is an array indicating elements which have a DAD defined
        DADisDefined = false(sizZ,sizX);
        
        % loop through bixels
        for zz = 1:sizZ
            for xx = 1:sizX
                
                % determine bixel coordinates
                x = X(xx);
                z = Z(zz);
                
                % insert coordinates into phase0 arrays
                DADx_phase0(zz,xx) = x;
                DADz_phase0(zz,xx) = z;
                
                % determine if there is a DAD defined at this bixels
                % eventually this should be compared against the identity
                % DAD, not DAD{1}
                [DADisDefined_temp,DADind] = ismember([x 0 z],masterRayPosBEV,'rows');
                
                if DADisDefined_temp
                    
                    % overwrite the temp variable if the DAD are NaNs
                    % (i.e., the DAD is not actually defined)
                    DADisDefined(zz,xx) = ~any(isnan(stf(i).DAD{phase}(DADind,:)));
                    
                    if DADisDefined(zz,xx)
                        % if the DAD is already defined, insert it into the
                        % array
                        stf(i).DADx{phase}(zz,xx) = stf(i).DAD{phase}(DADind,1);
                        stf(i).DADz{phase}(zz,xx) = stf(i).DAD{phase}(DADind,3);
                    end
                elseif xx == 1 || xx == sizX || zz == 1 || zz == sizZ
                    % if we are at the edge of the array, insert the
                    % identity
                    stf(i).DADx{phase}(zz,xx) = x;
                    stf(i).DADz{phase}(zz,xx) = z;
                    
                    % also, track that we have defined DAD here
                    DADisDefined(zz,xx) = true;
                end
            end
        end
        
        % now do interpolation to find the rest of the DAD values
        stf(i).DADx{phase}(~DADisDefined) = griddata(DADx_phase0(DADisDefined),DADz_phase0(DADisDefined),stf(i).DADx{phase}(DADisDefined),DADx_phase0(~DADisDefined),DADz_phase0(~DADisDefined));
        stf(i).DADz{phase}(~DADisDefined) = griddata(DADx_phase0(DADisDefined),DADz_phase0(DADisDefined),stf(i).DADz{phase}(DADisDefined),DADx_phase0(~DADisDefined),DADz_phase0(~DADisDefined));
        
    end
    
    matRad_progress(i,numel(stf));
end
