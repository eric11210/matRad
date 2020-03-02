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

if pln.propOpt.run4D && ~pln.propOpt.VMAToptions.fixedGantrySpeed
    warning('For 4D optimization, only fixed gantry speed is supported!');
    
    % set fixedGantrySpeed to true
    pln.propOpt.VMAToptions.fixedGantrySpeed = true;
end

% calculate total angular range
angularRange = abs(pln.propOpt.VMAToptions.finishingAngle-pln.propOpt.VMAToptions.startingAngle);

if pln.propOpt.VMAToptions.continuousAperture
    
    if pln.propOpt.run4D
        % angularRange = fluGantryAngleSpacing*numFluGantryAngles
        % also, deliveryTime = numFluGantryAngles*deltaT_sample
        numFluGantryAngles = round(pln.propOpt.VMAToptions.deliveryTime/pln.propOpt.prop4D.motionModel.deltaT_sample);
        fluGantryAngleSpacing = angularRange./numFluGantryAngles;
    else
        % angularRange = fluGantryAngleSpacing*numFluGantryAngles
        % ensure that fluGantryAngleSpacing < maxFluGantryAngleSpacing (as 
        % close as possible)
        numFluGantryAngles = ceil(angularRange./pln.propOpt.VMAToptions.maxFluGantryAngleSpacing);
        fluGantryAngleSpacing = angularRange./numFluGantryAngles;
    end
    
    % numGantryAngles*gantryAngleSpacing = numFluGantryAngles*fluGantryAngleSpacing
    % where
    % ensure that gantryAngleSpacing < maxgantryAngleSpacing (as close as
    % possible)
    numGantryAngles = ceil(numFluGantryAngles.*fluGantryAngleSpacing./pln.propOpt.VMAToptions.maxGantryAngleSpacing);
    % now ensure that numFluGantryAngles is an odd multiple of numGantryAngles so
    % that they align
    if pln.propOpt.run4D
        numFluGantryAngles = numGantryAngles.*(round(numFluGantryAngles./(2.*numGantryAngles)).*2+1);
    else
        numFluGantryAngles = numGantryAngles.*(floor(numFluGantryAngles./(2.*numGantryAngles)).*2+1);
    end
    gantryAngleSpacing = angularRange/numGantryAngles;
    
    
    % (numDAOGantryAngles-1)*DAOGantryAngleSpacing = (numGantryAngles-1)*gantryAngleSpacing
    % where
    % ensure that DAOGantryAngleSpacing < maxDAOGantryAngleSpacing (as close as
    % possible)
    numDAOGantryAngles = ceil((numGantryAngles-1).*gantryAngleSpacing./pln.propOpt.VMAToptions.maxDAOGantryAngleSpacing)+1;
    % now ensure that numGantryAngles-1 is a multiple of numDAOGantryAngles-1 so
    % that they align
    numGantryAngles = (numDAOGantryAngles-1).*ceil((numGantryAngles-1)./(numDAOGantryAngles-1))+1;
    %gantryAngleSpacing = angularRange./numGantryAngles;
    %DAOGantryAngleSpacing = (angularRange-gantryAngleSpacing)/(numDAOGantryAngles-1);
    % now ensure that numFluGantryAngles is an odd multiple of numGantryAngles so
    % that they align
    if pln.propOpt.run4D
        numFluGantryAngles = numGantryAngles.*(round(numFluGantryAngles./(2.*numGantryAngles)).*2+1);
    else
        numFluGantryAngles = numGantryAngles.*(floor(numFluGantryAngles./(2.*numGantryAngles)).*2+1);
    end
    fluGantryAngleSpacing = angularRange/numFluGantryAngles;
    gantryAngleSpacing = angularRange/numGantryAngles;
    DAOGantryAngleSpacing = (angularRange-gantryAngleSpacing)/(numDAOGantryAngles-1);
    
    % first and last gantry angles are in centre of arc
    firstGantryAngle = pln.propOpt.VMAToptions.startingAngle+gantryAngleSpacing/2;
    lastGantryAngle = pln.propOpt.VMAToptions.finishingAngle-gantryAngleSpacing/2;
    
    firstFluGantryAngle = pln.propOpt.VMAToptions.startingAngle+fluGantryAngleSpacing/2;
    lastFluGantryAngle = pln.propOpt.VMAToptions.finishingAngle-fluGantryAngleSpacing/2;
    
else
    
    % ensure that an integer number of DAO gantry angles will fit in angularRange degrees,
    % spaced at least as close as maxDAOGantryAngleSpacing
    numDAOGantryAngles = ceil(angularRange/pln.propOpt.VMAToptions.maxDAOGantryAngleSpacing);
    DAOGantryAngleSpacing = angularRange/numDAOGantryAngles;
    % actual number of gantry angles is numDAOGantryAngles+1;
    % ensure that DAOGantryAngleSpacing is an integer multiple of gantryAngleSpacing
    numGantryAngles = ceil(numDAOGantryAngles*DAOGantryAngleSpacing/pln.propOpt.VMAToptions.maxGantryAngleSpacing);
    gantryAngleSpacing = angularRange/numGantryAngles;
    
    % first and last gantry angles are at beginning and end of arc
    firstGantryAngle = pln.propOpt.VMAToptions.startingAngle;
    lastGantryAngle = pln.propOpt.VMAToptions.finishingAngle;
end

% calculate fluence angle spacing

% ensure that FMOGantryAngleSpacing is an odd integer multiple of DAOGantryAngleSpacing
numApertures = floor(pln.propOpt.VMAToptions.maxFMOGantryAngleSpacing/DAOGantryAngleSpacing);
if mod(numApertures,2) == 0
    numApertures = numApertures-1;
end
FMOGantryAngleSpacing = numApertures*DAOGantryAngleSpacing;

firstFMOGantryAngle = firstGantryAngle+DAOGantryAngleSpacing*floor(numApertures/2);
lastFMOGantryAngle = lastGantryAngle-DAOGantryAngleSpacing*floor(numApertures/2);


% insert final angle spacings
pln.propOpt.VMAToptions.fluGantryAngleSpacing   = fluGantryAngleSpacing;
pln.propOpt.VMAToptions.gantryAngleSpacing      = gantryAngleSpacing;
pln.propOpt.VMAToptions.DAOGantryAngleSpacing   = DAOGantryAngleSpacing;
pln.propOpt.VMAToptions.FMOGantryAngleSpacing   = FMOGantryAngleSpacing;

% define angles
pln.propStf.fluGantryAngles = firstFluGantryAngle:fluGantryAngleSpacing:lastFluGantryAngle;
pln.propStf.gantryAngles    = firstGantryAngle:gantryAngleSpacing:lastGantryAngle;
pln.propStf.DAOGantryAngles = firstGantryAngle:DAOGantryAngleSpacing:lastGantryAngle;
pln.propStf.FMOGantryAngles = firstFMOGantryAngle:FMOGantryAngleSpacing:lastFMOGantryAngle;

% everything else
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.couchAngles     = 0*pln.propStf.gantryAngles;
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

% update delivery time
if pln.propOpt.run4D
    pln.propOpt.VMAToptions.deliveryTime = pln.propOpt.prop4D.motionModel.deltaT_sample*numFluGantryAngles;
end



