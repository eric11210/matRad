function recalc = matRad_doseRecalc(cst,pln,recalc,ct,apertureInfo,calcDoseDirect,dij,doseOrFlu)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to recalculate dose on a Dij angular resolution, or using
% the dynamic method, whichever the user wants%
%
% call
%   resultGUI = matRad_doseRecalc(dij,apertureInfo,resultGUI,pln)
%
% input
%   dij:            matRad dij struct
%   apertureInfo:   aperture shape info struct
%   resultGUI:      resultGUI struct to which the output data will be added, if
%                   this field is empty optResult struct will be created
%                   (optional)
%
% output
%   resultGUI:  struct containing optimized fluence vector, dose, and
%               shape info
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

if nargin < 6
    calcDoseDirect = true;
end

if isfield(apertureInfo,'scaleFacRx')
    %weights were scaled to acheive 95% PTV coverage
    %scale back to "optimal" weights
    apertureInfo.apertureVector(1:apertureInfo.totalNumOfShapes*apertureInfo.numPhases) = apertureInfo.apertureVector(1:apertureInfo.totalNumOfShapes*apertureInfo.numPhases)/apertureInfo.scaleFacRx;
    
    % update aperture info vector
    apertureInfo = matRad_daoVec2ApertureInfo_bixWeightOnly(apertureInfo,apertureInfo.apertureVector);
end

recalc.apertureInfo = apertureInfo;

%recalculate dose with finer gantry angles
%to do this, we need Dij matrices at these new angles
%Calculate dose directly

%first, we need to update/generate new apertures for these angles
run4D                       = recalc.pln.propOpt.run4D;
recalc.pln.propOpt.run4D    = false;
recalc.pln                  = matRad_VMATGantryAngles(recalc.pln,cst,ct);
recalc.pln.propOpt.run4D    = run4D;
if ~recalc.interpNew
    %we will not interpolate new apertures
    %easiest way to to this is to make ALL gantryAngles optGantryAngles
    recalc.pln.propStf.DAOGantryAngles = recalc.pln.propStf.gantryAngles;
end

cd stf

switch doseOrFlu
    case 'dose'
        fname = sprintf('%.3f deg.mat',recalc.pln.propOpt.VMAToptions.maxGantryAngleSpacing);
    case 'flu'
        fname = sprintf('%.3f deg.mat',recalc.pln.propOpt.VMAToptions.maxFluGantryAngleSpacing);
end

if exist(fname,'file')
    load(fname,'stf');
else
    stf = matRad_generateStf(ct,cst,recalc.pln);
    save(fname,'stf')
end
recalc.stf = stf;
clear stf
cd ..

if ~recalc.interpNew || ~recalc.dijNew
    %duplicate any beam angles that are directly between two old
    %ones
    duplicate = false(size(recalc.pln.propStf.gantryAngles));
    for i = 1:numel(recalc.pln.propStf.gantryAngles)
        if numel(find(abs(recalc.pln.propStf.gantryAngles(i)-pln.propStf.gantryAngles) == min(abs(recalc.pln.propStf.gantryAngles(i)-pln.propStf.gantryAngles)))) > 1
            duplicate(i) = true;
        end
    end
    newGantryAngles = zeros(1,numel(recalc.pln.propStf.gantryAngles)+nnz(duplicate));
    newCouchAngles = zeros(1,numel(recalc.pln.propStf.gantryAngles)+nnz(duplicate));
    tempStf = recalc.stf;
    recalc.stf(1).copyInd = [];
    tempStf(1).copyInd = [];
    recalc.stf(1).stfCorr = [];
    tempStf(1).stfCorr = [];
    j = 1;
    for i = 1:numel(recalc.pln.propStf.gantryAngles)
        if duplicate(i)
            tempStf(j).stfCorr = false;
            newGantryAngles(j) = recalc.pln.propStf.gantryAngles(i);
            newCouchAngles(j) = recalc.pln.propStf.couchAngles(i);
            tempStf(j) = recalc.stf(i);
            tempStf(j).gantryAngle = recalc.stf(i-1).gantryAngle;
            tempStf(j).copyInd = 1;
            
            j = j+1;
            
            newGantryAngles(j) = recalc.pln.propStf.gantryAngles(i);
            newCouchAngles(j) = recalc.pln.propStf.couchAngles(i);
            tempStf(j) = recalc.stf(i);
            tempStf(j).gantryAngle = recalc.stf(i+1).gantryAngle;
            tempStf(j).copyInd = 2;
        else
            tempStf(j).stfCorr = true;
            newGantryAngles(j) = recalc.pln.propStf.gantryAngles(i);
            newCouchAngles(j) = recalc.pln.propStf.couchAngles(i);
            tempStf(j) = recalc.stf(i);
        end
        j = j+1;
    end
    recalc.pln.propStf.gantryAngles = newGantryAngles;
    recalc.pln.propStf.couchAngles = newCouchAngles;
    recalc.pln.propStf.numOfBeams = numel(recalc.pln.propStf.gantryAngles);
    %recalc.pln.optGantryAngles = recalc.pln.gantryAngles;
    recalc.stf = tempStf;
end

if pln.propOpt.VMAToptions.continuousAperture
    recalc = matRad_recalcApertureInfoFluAngles(recalc,recalc.apertureInfo);
else
    recalc = matRad_recalcApertureInfoFromStatic(recalc,recalc.apertureInfo);
end

if ~recalc.interpNew || ~recalc.dijNew
    tempPln = recalc.pln;
    tempStf = recalc.stf;
    for i = 1:numel(tempPln.propStf.gantryAngles)
        diff = abs(tempPln.propStf.gantryAngles(i)-pln.propStf.gantryAngles);
        minDiffInd = diff == min(diff);
        minDiffInd1 = find(tempPln.propStf.gantryAngles == pln.propStf.gantryAngles(find(minDiffInd,1,'first')));
        minDiffInd2 = find(tempPln.propStf.gantryAngles == pln.propStf.gantryAngles(find(minDiffInd,1,'last')));
        
        if ~recalc.dijNew
            if isempty(recalc.stf(i).copyInd)
                recalc.stf(i) = tempStf(minDiffInd1);
                recalc.pln.propStf.gantryAngles(i) = tempPln.propStf.gantryAngles(minDiffInd1);
            elseif recalc.stf(i).copyInd == 1
                recalc.stf(i) = tempStf(minDiffInd1);
                recalc.pln.propStf.gantryAngles(i) = tempPln.propStf.gantryAngles(minDiffInd1);
            elseif recalc.stf(i).copyInd == 2
                recalc.stf(i) = tempStf(minDiffInd2);
                recalc.pln.propStf.gantryAngles(i) = tempPln.propStf.gantryAngles(minDiffInd2);
            end
        elseif ~recalc.interpNew
            if numel(minDiffInd) > 1
                recalc.stf(i).gantryAngle = tempPln.propStf.gantryAngles(i);
            end
        end
    end
end

if calcDoseDirect
    clear global
    recalc.resultGUI = matRad_calcDoseDirect(ct,recalc.stf,recalc.pln,cst,recalc.apertureInfo.bixelWeights);
else
    recalc.resultGUI.w = recalc.apertureInfo.bixelWeights;
    
    % set dose calculation options
    options.bioOpt          = recalc.pln.propOpt.bioOptimization;
    if recalc.pln.propOpt.run4D
        options.numOfScenarios = dij.numPhases;
    else
        options.numOfScenarios  = dij.numOfScenarios;
    end
    dij.scaleFactor = recalc.apertureInfo.weightToMU./dij.weightToMU;
    d = matRad_backProjection(recalc.resultGUI.w,dij,options);
    recalc.resultGUI.physicalDose = reshape(d,dij.dimensions);
    
    % adjust overlap priorities
    cst_Over = matRad_setOverlapPriorities(cst);
    
    % adjust objectives _and_ constraints internally for fractionation
    for i = 1:size(cst_Over,1)
        for j = 1:size(cst_Over{i,6},1)
            cst_Over{i,6}(j).dose = cst_Over{i,6}(j).dose/recalc.pln.numOfFractions;
        end
    end
    
    % bixel based objective function calculation
    recalc.f = matRad_objFuncWrapper(recalc.resultGUI.w,dij,cst_Over,options);
end
