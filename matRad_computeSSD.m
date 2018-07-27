function stf = matRad_computeSSD(stf,ct,mode)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad SSD calculation
% 
% call
%   stf = matRad_computeSSD(stf,ct,mode)
%
% input
%   ct:             ct cube
%   stf:            matRad steering information struct
%   mode:           optional parameter specifying how to handle multiple
%                   cubes to compute one SSD
% output
%   stf:            matRad steering information struct
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    mode = 'first';
end

% booleon to show warnings only once in the console
boolShowWarning = true;

% set density threshold for SSD computation
densityThreshold = 0.05;

if strcmp(mode,'first')
    
    fprintf('matRad: SSD calculation... ');
    
    for i = 1:size(stf,2)
        SSD = cell(ct.numOfCtScen,stf(i).numOfRays);
        for j = 1:stf(i).numOfRays
            [alpha,~,rho,~,~] = matRad_siddonRayTracer(stf(i).isoCenter, ...
                                 ct.resolution, ...
                                 stf(i).sourcePoint, ...
                                 stf(i).ray(j).targetPoint, ...
                                 ct.cube);
                             
            for k = 1:ct.numOfCtScen
                ixSSD = find(rho{k} > densityThreshold,1,'first');
                
                if boolShowWarning
                    if isempty(ixSSD)
                        warning('ray does not hit patient. Trying to fix afterwards...');
                        boolShowWarning = false;
                    elseif ixSSD(1) == 1
                        warning('Surface for SSD calculation starts directly in first voxel of CT\n');
                        boolShowWarning = false;
                    end
                end
                
                % calculate SSD
                SSD{k,j} = double(2 * stf(i).SAD * alpha(ixSSD));
                stf(i).ray(j).SSD{k} = SSD{k,j};
            end
        end
        
        % try to fix SSD by using SSD of closest neighbouring ray
        SSDnotSet = find(cellfun('isempty',SSD))';
        if ~isempty(SSDnotSet)
            rayPos_bev = reshape([stf(i).ray(:).rayPos_bev]',[3 stf(i).numOfRays])';
            for ind = SSDnotSet
                [k, j] = ind2sub(size(SSD),ind);
                stf(i).ray(j).SSD{k} =  matRad_closestNeighbourSSD(rayPos_bev, SSD(k,:), rayPos_bev(j,:));
            end
        end
        
        matRad_progress(i,size(stf,2));
    end
else
    error('mode not defined for SSD calculation');
end


% default setting only use first cube
function bestSSD = matRad_closestNeighbourSSD(rayPos, SSD, currPos)
    vDistances = sum((rayPos - repmat(currPos,size(rayPos,1),1)).^2,2);
    [~, vIdx]   = sort(vDistances);
    for ix = vIdx'
        bestSSD = SSD{ix};
        % if SSD has been found, bestSSD is not empty
        if ~any(isempty(bestSSD))
            break
        end
    end
    if any(isempty(bestSSD))
        error('Could not fix SSD calculation.');
    end
end






end
