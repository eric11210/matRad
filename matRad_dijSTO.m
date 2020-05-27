function [dij_STO,trajectory] = matRad_dijSTO(dij,pln,stf,apertureInfo)

%% select most probable trajectory

% prepare model
motionModel = matRad_prepModelForOpt(pln,stf,apertureInfo);

% determine transition times
tTrans = [apertureInfo.beam.time]';

% determine the times for each step, which are the same for all histories
tSimulated = [0; cumsum(tTrans)];

% calculate the number of steps to take (it's possible that there is really
% one more step than necessary, but it's good to have more steps than we
% require than not enough)
nSteps = 1+ceil(tSimulated(end)./motionModel.deltaT_sample);

% sample the most probable trajectory
lSimulated = matRad_runMarkovChain_P(motionModel,nSteps,true);

% convert subphase to phase
pSimulated = motionModel.indices.subPhase2PosPhase(lSimulated);

% store trajectory information
trajectory.tSimulated = tSimulated;
trajectory.pSimulated = pSimulated;

%% construct effective dij

dij_fieldNames = fieldnames(dij);

for i = 1:numel(dij_fieldNames)
    
    if ~strcmp(dij_fieldNames{i},'physicalDose')
        
        dij_STO.(dij_fieldNames{i}) = dij.(dij_fieldNames{i});
    end
end

dij_STO.physicalDose{1} = spalloc(dij.numOfVoxels,dij.totalNumOfBixels,1);
dij.numOfScenarios      = 1;


for i = 1:numel(stf)
    
    % get last and next bixel index maps
    lastBixelIndMap = apertureInfo.beam(i).lastBixelIndMap(~isnan(apertureInfo.beam(i).lastBixelIndMap));
    nextBixelIndMap = apertureInfo.beam(i).nextBixelIndMap(~isnan(apertureInfo.beam(i).nextBixelIndMap));
    
    % give dose to last beam
    % initial phase
    dij_STO.physicalDose{1}(:,lastBixelIndMap) = dij_STO.physicalDose{1}(:,lastBixelIndMap) + ...
        stf(i).propVMAT.absFracToLastDose_arcI.*dij.physicalDose{pSimulated(i)}(:,lastBixelIndMap);
    % final phase
    dij_STO.physicalDose{1}(:,lastBixelIndMap) = dij_STO.physicalDose{1}(:,lastBixelIndMap) + ...
        stf(i).propVMAT.absFracToLastDose_arcF.*dij.physicalDose{pSimulated(i+1)}(:,lastBixelIndMap);
    
    % give dose to next beam
    % initial phase
    dij_STO.physicalDose{1}(:,nextBixelIndMap) = dij_STO.physicalDose{1}(:,nextBixelIndMap) + ...
        stf(i).propVMAT.absFracToNextDose_arcI.*dij.physicalDose{pSimulated(i)}(:,nextBixelIndMap);
    % final phase
    dij_STO.physicalDose{1}(:,nextBixelIndMap) = dij_STO.physicalDose{1}(:,nextBixelIndMap) + ...
        stf(i).propVMAT.absFracToNextDose_arcF.*dij.physicalDose{pSimulated(i+1)}(:,nextBixelIndMap);
    
end

end