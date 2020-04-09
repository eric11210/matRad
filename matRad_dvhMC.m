function [pdvh_MC,dvh_mean_MC,dvh_std_MC,dMean_MC,dVar_MC] = matRad_dvhMC(apertureInfo,dij,cst,pln,nHistories,p)

% shuffle rng
rng('shuffle');

% initialize options
if apertureInfo.run4D
    options.numOfScenarios = dij.numPhases;
else
    options.numOfScenarios  = dij.numOfScenarios;
end
options.bioOpt = 'none';
options.run4D = apertureInfo.run4D;
options.FMO = false;

% default percentiles
if nargin < 6
    p = [5 25 50 75 95];
end

% initialize sum of dose and square of dose
dSum    = zeros(dij.numOfVoxels,1);
d2Sum   = zeros(dij.numOfVoxels,1);

% get dose for mean
d_mean = matRad_backProjection(apertureInfo.bixelWeights,dij,options);

% now the number of scenarios should be the actual number of phases
options.numOfScenarios = dij.numPhases;

% reshape dose
dCubes_Mean = reshape(d_mean,dij.dimensions);

% get DVH mean
dvh_mean = matRad_calcDVH(cst,dCubes_Mean,'cum');

dvh_mean_unique = dvh_mean;
dvh_sum         = dvh_mean;
dvh_2sum        = dvh_mean;
dvh_hists       = dvh_mean;

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
    
    % prepare DVH for all histories
    dvh_hists(struct).volumePoints  = zeros(nHistories,1000);
    dvh_hists(struct).doseGrid      = zeros(nHistories,1000);
    
end

% determine transition times
tTrans = [apertureInfo.beam.time]';

% prepare motionModel structure
apertureInfo_frac                       = apertureInfo;
apertureInfo_frac                       = rmfield(apertureInfo_frac,'motionModel');
apertureInfo_frac.motionModel.type      = 'single';
apertureInfo_frac.motionModel.numPhases = apertureInfo_frac.numPhases;
apertureInfo_frac.motionModel.indices.nSubPhases = apertureInfo_frac.numPhases;

[apertureInfo_frac.motionModel.indices.subPhase2PosPhase_gridJ,apertureInfo_frac.motionModel.indices.subPhase2PosPhase_gridI] = meshgrid(1:apertureInfo_frac.numPhases);

% determine the times for each step, which are the same for all histories
tSimulated_hist = [0; cumsum(tTrans)];

% calculate the number of steps to take (it's possible that there is really
% one more step than necessary, but it's good to have more steps than we
% require than not enough)
nSteps = 1+ceil(tSimulated_hist(end)./apertureInfo.motionModel.deltaT_sample);

%% do MC

% loop over number of histories
for hist = 1:nHistories
    
    % initialize dose for history (sum over fractions)
    d_hist = zeros(dij.numOfVoxels,1);
    
    % loop over number of fractions
    for frac = 1:pln.numOfFractions
        
        % MC simulation of tumour trajectory
        lSimulated_deltaT_hist = matRad_runMarkovChain_P(apertureInfo.motionModel,nSteps);
        
        % interpolate subphase trajectory to desired temporal frequency, convert
        % to posPhase
        lSimulated_hist = matRad_resampleTrajAndSubPhase2PosPhase(lSimulated_deltaT_hist,tSimulated_hist,apertureInfo.motionModel);
        
        % insert trajectory in apertureInfo struct
        apertureInfo_frac.motionModel.lSimulated = lSimulated_hist;
        apertureInfo_frac.motionModel.tSimulated = tSimulated_hist;
        
        % get weights for this particular history
        apertureInfo_frac = matRad_daoVec2ApertureInfo_bixWeightOnly(apertureInfo_frac,apertureInfo_frac.apertureVector);
        
        % get dose for this particular trajectory
        d_frac = matRad_backProjection(apertureInfo_frac.bixelWeights,dij,options);
        
        % sum fraction dose into total dose
        d_hist = d_hist+d_frac;
        
        % update dSum and d2Sum
        dSum    = dSum + d_frac;
        d2Sum   = d2Sum + d_frac.^2;
    end
    
    % reshape dose
    dCubes_hist = reshape(d_hist,dij.dimensions);
    
    % get DVH for history
    dvh_hist = matRad_calcDVH(cst,dCubes_hist,'cum');
    
    for struct = 1:numel(dvh_hist)
        
        [volumePointsUnique,iOrig,~]    = unique(dvh_hist(struct).volumePoints,'stable');
        doseGridUnique                  = dvh_hist(struct).doseGrid(iOrig);
        
        doseGrid_hist = interp1(volumePointsUnique,doseGridUnique,dvh_mean_unique(struct).volumePoints);
        
        %doseGrid_hist = interp1(dvh_mean_unique(struct).volumePoints,dvh_mean_unique(struct).doseGrid,dvh_hist(struct).volumePoints);
        
        
        % calculate sum of dose and sum of square dose
        dvh_sum(struct).doseGrid    = dvh_sum(struct).doseGrid+doseGrid_hist;
        dvh_2sum(struct).doseGrid   = dvh_2sum(struct).doseGrid+doseGrid_hist.^2;
        
        % insert DVH for history into DVH for all histories
        dvh_hists(struct).doseGrid(hist,:)      = dvh_hist(struct).doseGrid;
        dvh_hists(struct).volumePoints(hist,:)  = dvh_hist(struct).volumePoints;
        
    end
    
end

% calculate mean and variance of dose
dMean_MC   = dSum./(nHistories.*pln.numOfFractions);
dVar_MC    = (d2Sum-nHistories.*pln.numOfFractions.*dMean_MC.^2)./(nHistories.*pln.numOfFractions-1);

%% calculate mean and std DVHs, and PDVHs

dvh_mean_MC = dvh_mean_unique;
dvh_std_MC  = dvh_mean_unique;

pdvh_MC = dvh_mean;

% construct percentile DVHs for each structure
for struct = 1:numel(dvh_mean)
    
    % get max dose for structure
    maxDose = max(dvh_hists(struct).doseGrid(:,end));
    
    % construct dose grid for structure
    pdvh_MC(struct).doseGrid = linspace(0,maxDose,1000);
    
    for hist = 1:nHistories
        % reinterpolate the DVH of each history to match this grid
        dvh_hists(struct).volumePoints(hist,:)  = interp1(dvh_hists(struct).doseGrid(hist,:),dvh_hists(struct).volumePoints(hist,:),pdvh_MC(struct).doseGrid);
        dvh_hists(struct).doseGrid(hist,:)      = pdvh_MC(struct).doseGrid;
        
    end
    
    % get the percentiles
    pdvh_MC(struct).volumePoints = prctile(dvh_hists(struct).volumePoints,p);
    
    dvh_mean_MC(struct).volumePoints    = dvh_mean_unique(struct).volumePoints;
    dvh_mean_MC(struct).doseGrid        = dvh_sum(struct).doseGrid./nHistories;
    dvh_mean_MC(struct).name            = dvh_mean(struct).name;
    
    dvh_std_MC(struct).volumePoints    = dvh_mean_unique(struct).volumePoints;
    dvh_std_MC(struct).doseGrid        = sqrt((dvh_2sum(struct).doseGrid-nHistories.*dvh_mean_MC(struct).doseGrid.^2)./(nHistories-1));
    %dvh_std_MC(struct).doseGrid        = sqrt((dvh_2sum(struct).doseGrid-dvh_mean_unique(struct).doseGrid.^2)./(nHistories));
    dvh_std_MC(struct).name            = dvh_mean(struct).name;
    
end

end

