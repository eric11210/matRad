function recalc = matRad_recalcApertureInfoFluAngles(recalc,apertureInfoOld)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to apertures for a different dose resolution
%
% call
%   recalc = matRad_recalcApertureInfo(recalc,apertureInfo)
%
% input
%   recalc:
%   apertureInfo:
%
% output
%   recalc:
%
% References
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

% strategy for this is to put in the minimum set of metadata into the
% apertureInfoNew structure, then use daoVec2ApertureInfo to get the leaf
% positions, bixel weights, etc. at the new angles

% extract current pln and stf structures
pln = recalc.pln;
stf = recalc.stf;

% save global data
apertureInfoNew.runVMAT             = pln.propOpt.runVMAT;
apertureInfoNew.preconditioner      = pln.propOpt.preconditioner;
apertureInfoNew.run4D               = pln.propOpt.run4D;
apertureInfoNew.varOpt              = pln.propOpt.varOpt;
apertureInfoNew.numPhases           = apertureInfoOld.numPhases;
apertureInfoNew.bixelWidth          = pln.propStf.bixelWidth;
apertureInfoNew.numOfMLCLeafPairs   = apertureInfoOld.numOfMLCLeafPairs;
apertureInfoNew.totalNumOfBixels    = sum([stf(:).totalNumOfBixels]);
apertureInfoNew.totalNumOfOptBixels = apertureInfoOld.totalNumOfOptBixels;
apertureInfoNew.totalNumOfShapes    = apertureInfoOld.totalNumOfShapes;
apertureInfoNew.weightToMU          = apertureInfoOld.weightToMU;

% preallocate propVMAT.beam (metadata)
%stf(stf(1).propVMAT.beamParentIndex).propVMAT.leafDir = 1;
apertureInfoNew.propVMAT.beam(size(stf,2)) = stf(stf(1).propVMAT.beamParentIndex).propVMAT;

% preallocate jacobT
apertureInfoNew.propVMAT.jacobT = zeros(apertureInfoNew.totalNumOfShapes,size(stf,2));

% get apertureVector
apertureVector = apertureInfoOld.apertureVector;

% loop over all beams
for i = 1:size(stf,2)
    
    % put all propVMAT stuff from stf into apertureInfo
    apertureInfoNew.propVMAT.beam(i) = stf(i).propVMAT;
    
    if stf(i).propVMAT.doseBeam
        
        % find the beam index from the old apertureInfo structure
        iOld = find(abs(stf(i).gantryAngle - [apertureInfoOld.beam.gantryAngle]) < 1e-8);
        
        % for fluence-only beams, we will determine bixel index maps later
        apertureInfoNew.beam(i).bixelIndMap = apertureInfoOld.beam(iOld).bixelIndMap;
        apertureInfoNew.beam(i).numBixels   = apertureInfoOld.beam(iOld).numBixels;
    end
    
    if stf(i).propVMAT.DAOBeam
        
        % loop over all phases
        for phase = 1:apertureInfoNew.numPhases
            
            % get new offset, jacobiScales
            apertureInfoNew.beam(i).shape{phase}.weightOffset   = apertureInfoOld.beam(iOld).shape{phase}.weightOffset;
            apertureInfoNew.beam(i).shape{phase}.vectorOffset   = apertureInfoOld.beam(iOld).shape{phase}.vectorOffset;
            apertureInfoNew.beam(i).shape{phase}.jacobiScale    = apertureInfoOld.beam(iOld).shape{phase}.jacobiScale;
            
            % rescale the weight to account for the change in fluAngleBordersDiff
            apertureInfoNew.beam(i).shape{phase}.weight = apertureInfoOld.beam(iOld).shape{phase}.weight.*apertureInfoNew.propVMAT.beam(i).fluAngleBordersDiff./apertureInfoOld.propVMAT.beam(iOld).fluAngleBordersDiff;
            
            % input new weight in apertureVector
            apertureVector(apertureInfoNew.beam(i).shape{phase}.weightOffset) = apertureInfoNew.beam(i).shape{phase}.weight.*apertureInfoNew.beam(i).shape{phase}.jacobiScale;
        end
        
        apertureInfoNew.propVMAT.jacobT(stf(i).propVMAT.DAOIndex,i) = stf(i).propVMAT.timeFacCurr;
    else
        apertureInfoNew.propVMAT.jacobT(stf(stf(i).propVMAT.lastDAOIndex).propVMAT.DAOIndex,i) = apertureInfoNew.propVMAT.jacobT(stf(stf(i).propVMAT.lastDAOIndex).propVMAT.DAOIndex,i)+stf(stf(i).propVMAT.lastDAOIndex).propVMAT.timeFacCurr.*stf(i).propVMAT.fracFromLastDAO_gantryRot.*stf(i).propVMAT.fluAngleBordersDiff./stf(stf(i).propVMAT.lastDAOIndex).propVMAT.fluAngleBordersDiff;
        apertureInfoNew.propVMAT.jacobT(stf(stf(i).propVMAT.nextDAOIndex).propVMAT.DAOIndex,i) = apertureInfoNew.propVMAT.jacobT(stf(stf(i).propVMAT.nextDAOIndex).propVMAT.DAOIndex,i)+stf(stf(i).propVMAT.nextDAOIndex).propVMAT.timeFacCurr.*stf(i).propVMAT.fracFromNextDAO_gantryRot.*stf(i).propVMAT.fluAngleBordersDiff./stf(stf(i).propVMAT.nextDAOIndex).propVMAT.fluAngleBordersDiff;
    end
    
    
    % save data for each beam
    apertureInfoNew.beam(i).numOfShapes = 1;
    apertureInfoNew.beam(i).numOfActiveLeafPairs = apertureInfoOld.beam(1).numOfActiveLeafPairs;
    apertureInfoNew.beam(i).leafPairPos = apertureInfoOld.beam(1).leafPairPos;
    apertureInfoNew.beam(i).isActiveLeafPair = apertureInfoOld.beam(1).isActiveLeafPair;
    apertureInfoNew.beam(i).centralLeafPair = apertureInfoOld.beam(1).centralLeafPair;
    apertureInfoNew.beam(i).lim_l = apertureInfoOld.beam(1).lim_l;
    apertureInfoNew.beam(i).lim_r = apertureInfoOld.beam(1).lim_r;
    apertureInfoNew.beam(i).posOfCornerBixel = apertureInfoOld.beam(1).posOfCornerBixel;
    apertureInfoNew.beam(i).MLCWindow = apertureInfoOld.beam(1).MLCWindow;
    apertureInfoNew.beam(i).gantryAngle = stf(i).gantryAngle;
end

% save more metadata
apertureInfoNew = matRad_apertureInfoMeta(apertureInfoNew,pln,stf,0);

% get updated leaf positions, bixel weights, etc. from apertureVector
apertureInfoNew = matRad_daoVec2ApertureInfo_bixWeightOnly(apertureInfoNew,apertureVector);

recalc.apertureInfo = apertureInfoNew;