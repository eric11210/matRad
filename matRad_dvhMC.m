function [dvh_mean_MC,dvh_std_MC] = matRad_dvhMC(apertureInfo,dij,cst,nHistories)

% initialize options
if apertureInfo.run4D
    options.numOfScenarios = dij.numPhases;
else
    options.numOfScenarios  = dij.numOfScenarios;
end
options.bioOpt = 'none';
options.run4D = apertureInfo.run4D;
options.FMO = false;

% get dose for mean
d_mean = matRad_backProjection(apertureInfo.bixelWeights,dij,options);

% reshape dose
dCubes_Mean = reshape(d_mean,dij.dimensions);

% get DVH mean
dvh_mean = matRad_calcDVH(cst,dCubes_Mean,'cum');

dvh_mean_unique = dvh_mean;
dvh_sum         = dvh_mean;
dvh_2sum        = dvh_mean;

for struct = 1:numel(dvh_mean)
    
    [volumePointsUnique,iOrig,~] = unique(dvh_mean(struct).volumePoints,'stable');
    
    dvh_mean_unique(struct).volumePoints    = volumePointsUnique;
    dvh_mean_unique(struct).doseGrid        = dvh_mean(struct).doseGrid(iOrig);
    dvh_mean_unique(struct).name            = dvh_mean(struct).name;
    
    dvh_sum(struct).volumePoints            = volumePointsUnique;
    dvh_sum(struct).doseGrid                = zeros(1,numel(volumePointsUnique));
    dvh_sum(struct).name                    = dvh_mean(struct).name;
    
    dvh_2sum(struct).volumePoints           = volumePointsUnique;
    dvh_2sum(struct).doseGrid               = zeros(1,numel(volumePointsUnique));
    dvh_2sum(struct).name                   = dvh_mean(struct).name;
    
end

% determine transition times
tTrans = [apertureInfo.beam.time]';

% prepare motionModel structure
apertureInfo_hist                       = apertureInfo;
apertureInfo_hist                       = rmfield(apertureInfo_hist,'motionModel');
apertureInfo_hist.motionModel.type      = 'single';
apertureInfo_hist.motionModel.numPhases = apertureInfo_hist.numPhases;
apertureInfo_hist.motionModel.indices.nSubPhases = apertureInfo_hist.numPhases;

[apertureInfo_hist.motionModel.indices.subPhase2PosPhase_gridJ,apertureInfo_hist.motionModel.indices.subPhase2PosPhase_gridI] = meshgrid(1:apertureInfo_hist.numPhases);

% loop over number of histories
for hist = 1:nHistories
    
    % MC simulation of tumour trajectory
    [lSimulated_hist,tSimulated_hist]= matRad_runMarkovChain_Q(apertureInfo.motionModel.qij,apertureInfo.motionModel.initProb,tTrans);
    
    % insert trajectory in apertureInfo struct
    apertureInfo_hist.motionModel.lSimulated = apertureInfo.motionModel.indices.subPhase2Phase(lSimulated_hist);
    apertureInfo_hist.motionModel.tSimulated = tSimulated_hist;
    
    % get weights for this particular history
    apertureInfo_hist = matRad_daoVec2ApertureInfo(apertureInfo_hist,apertureInfo_hist.apertureVector);
    
    % get dose for this particular trajectory
    d_hist = matRad_backProjection(apertureInfo_hist.bixelWeights,dij,options);
    
    % reshape dose
    dCubes_hist = reshape(d_hist,dij.dimensions);
    
    % get DVH for history
    dvh_hist = matRad_calcDVH(cst,dCubes_hist,'cum');
    
    for struct = 1:numel(dvh_hist)
        
        [volumePointsUnique,iOrig,~]    = unique(dvh_hist(struct).volumePoints,'stable');
        doseGridUnique                  = dvh_hist(struct).doseGrid(iOrig);
        
        doseGrid_hist = interp1(volumePointsUnique,doseGridUnique,dvh_mean_unique(struct).volumePoints);
        
        %doseGrid_hist = interp1(dvh_mean_unique(struct).volumePoints,dvh_mean_unique(struct).doseGrid,dvh_hist(struct).volumePoints);
        
        dvh_sum(struct).doseGrid    = dvh_sum(struct).doseGrid+doseGrid_hist;
        dvh_2sum(struct).doseGrid   = dvh_2sum(struct).doseGrid+doseGrid_hist.^2;
        
    end
end

dvh_mean_MC = dvh_mean_unique;
dvh_std_MC  = dvh_mean_unique;

for struct = 1:numel(dvh_mean)
    
    dvh_mean_MC(struct).volumePoints    = dvh_mean_unique(struct).volumePoints;
    dvh_mean_MC(struct).doseGrid        = dvh_sum(struct).doseGrid./nHistories;
    dvh_mean_MC(struct).name            = dvh_mean(struct).name;
    
    dvh_std_MC(struct).volumePoints    = dvh_mean_unique(struct).volumePoints;
    dvh_std_MC(struct).doseGrid        = sqrt((dvh_2sum(struct).doseGrid-nHistories.*dvh_mean_MC(struct).doseGrid.^2)./(nHistories-1));
    %dvh_std_MC(struct).doseGrid        = sqrt((dvh_2sum(struct).doseGrid-dvh_mean_unique(struct).doseGrid.^2)./(nHistories));
    dvh_std_MC(struct).name            = dvh_mean(struct).name;
    
end

end

