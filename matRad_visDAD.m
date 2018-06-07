function stf = matRad_visDAD(stf)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visuliaze DAD info in stf
%
% call
%   stf = matRad_visDAD(stf)
%
% input
%   ct:         ct cube
%   cst:        matRad cst struct
%   pln:        matRad plan meta information struct
%   visMode:    toggle on/off different visualizations by setting this value to 1,2,3 (optional)
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

for i = 1:numel(stf)
    
    DAD_phase1 = stf(i).DAD{1};
    figure
    
    for phase = 2:numel(stf(i).DAD)
        
        vel_phaseP = stf(i).DAD{phase}-stf(i).DAD{1};
        
        subplot(2,numel(stf(i).DAD)/2,phase)
        quiver(DAD_phase1(:,1),DAD_phase1(:,3),vel_phaseP(:,1),vel_phaseP(:,3));
        title(sprintf('Phase %d',phase));
        xlabel('horiz. pos. [mm]');
        ylabel('vert. pos. [mm]');
        
    end
    suptitle(sprintf('Gantry angle: %.1f%c',stf(i).gantryAngle,char(176)));
    
end