function pln = matRad_VMATGantryAngles(pln,cst,ct)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad determine gantry angles for VMAT
% 
% call
%   matRad_VMATGantryAngles(pln)
%
% input
%   pln:                matRad plan meta information struct
%
% output
%   pln:                matRad plan meta information struct
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

%should have already defined fields:
%   maxGantryAngleSpacing
%   maxDAOGantryAngleSpacing
%   maxFMOGantryAngleSpacing


if ~isfield(pln.propOpt.VMAToptions,'maxGantryAngleSpacing')
    error('Please define pln.propOpt.maxGantryAngleSpacing.');
end

if ~isfield(pln.propOpt.VMAToptions,'maxDAOGantryAngleSpacing')
    error('Please define pln.propOpt.maxDAOGantryAngleSpacing.');
end

if ~isfield(pln.propOpt.VMAToptions,'maxFMOGantryAngleSpacing')
    error('Please define pln.propOpt.maxFMOGantryAngleSpacing.');
end

angularRange = abs(pln.propOpt.VMAToptions.finishingAngle-pln.propOpt.VMAToptions.startingAngle);

% angularRange = gantryAngleSpacing*numGantryAngles
% ensure that gantryAngleSpacing < maxGantryAngleSpacing (as close as
% possible)
numGantryAngles = ceil(angularRange./pln.propOpt.VMAToptions.maxGantryAngleSpacing)-1;
gantryAngleSpacing = angularRange./(numGantryAngles+1);

% numDAOGantryAngles*DAOGantryAngleSpacing = numGantryAngles*gantryAngleSpacing
% where 
% ensure that DAOGantryAngleSpacing < maxDAOGantryAngleSpacing (as close as
% possible)
numDAOGantryAngles = ceil(numGantryAngles.*gantryAngleSpacing./pln.propOpt.VMAToptions.maxDAOGantryAngleSpacing);
% now ensure that numGantryAngles is a multiple of numDAOGantryAngles so
% that they align
numGantryAngles = numDAOGantryAngles.*ceil(numGantryAngles./numDAOGantryAngles);
gantryAngleSpacing = angularRange./(numGantryAngles+1);
DAOGantryAngleSpacing = (angularRange-gantryAngleSpacing)/numDAOGantryAngles;

firstGantryAngle = pln.propOpt.VMAToptions.startingAngle+gantryAngleSpacing/2;
lastGantryAngle = pln.propOpt.VMAToptions.finishingAngle-gantryAngleSpacing/2;

% ensure that FMOGantryAngleSpacing is an odd integer multiple of DAOGantryAngleSpacing
numApertures = floor(pln.propOpt.VMAToptions.maxFMOGantryAngleSpacing/DAOGantryAngleSpacing);
if mod(numApertures,2) == 0
    numApertures = numApertures-1;
end
FMOGantryAngleSpacing = numApertures*DAOGantryAngleSpacing;

firstFMOGantryAngle = firstGantryAngle+DAOGantryAngleSpacing*floor(numApertures/2);
lastFMOGantryAngle = lastGantryAngle-DAOGantryAngleSpacing*floor(numApertures/2);

% define angles
pln.propStf.gantryAngles    = firstGantryAngle:gantryAngleSpacing:lastGantryAngle;
pln.propStf.DAOGantryAngles = firstGantryAngle:DAOGantryAngleSpacing:lastGantryAngle;
pln.propStf.FMOGantryAngles = firstFMOGantryAngle:FMOGantryAngleSpacing:lastFMOGantryAngle;

% everything else
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.couchAngles     = 0*pln.propStf.gantryAngles;
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);



