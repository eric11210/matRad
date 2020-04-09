function [dMean_MC,dVar_MC] = matRad_doseVarianceMC(apertureInfo,dij,nHistories)
% Monte Carlo estimate of the mean and variance of the dose

% initialize sum of dose and square of dose
dSum    = zeros(dij.numOfVoxels,1);
d2Sum   = zeros(dij.numOfVoxels,1);

% determine transition times
tTrans = [apertureInfo.beam.time]';

% initialize options
if apertureInfo.run4D
    options.numOfScenarios = dij.numPhases;
else
    options.numOfScenarios  = dij.numOfScenarios;
end
options.bioOpt = 'none';
options.run4D = apertureInfo.run4D;
options.FMO = false;

% prepare motionModel structure
apertureInfo_hist                       = apertureInfo;
apertureInfo_hist                       = rmfield(apertureInfo_hist,'motionModel');
apertureInfo_hist.motionModel.type      = 'single';
apertureInfo_hist.motionModel.numPhases = apertureInfo_hist.numPhases;
apertureInfo_hist.motionModel.indices.nSubPhases = apertureInfo_hist.numPhases;

[apertureInfo_hist.motionModel.indices.subPhase2PosPhase_gridJ,apertureInfo_hist.motionModel.indices.subPhase2PosPhase_gridI] = meshgrid(1:apertureInfo_hist.numPhases);

% determine the times for each step, which are the same for all histories
tSimulated_hist = [0; cumsum(tTrans)];

% calculate the number of steps to take (it's possible that there is really
% one more step than necessary, but it's good to have more steps than we
% require than not enough)
nSteps = 1+ceil(tSimulated_hist(end)./apertureInfo.motionModel.deltaT_sample);

% loop over number of histories
for hist = 1:nHistories
    
    % MC simulation of tumour trajectory
    lSimulated_deltaT_hist = matRad_runMarkovChain_P(apertureInfo.motionModel,nSteps);
    
    % interpolate subphase trajectory to desired temporal frequency, convert
    % to posPhase
    lSimulated_hist = matRad_resampleTrajAndSubPhase2PosPhase(lSimulated_deltaT_hist,tSimulated_hist,apertureInfo.motionModel);
    
    % insert trajectory in apertureInfo struct
    apertureInfo_hist.motionModel.lSimulated = lSimulated_hist;
    apertureInfo_hist.motionModel.tSimulated = tSimulated_hist;
    
    % get weights for this particular history
    apertureInfo_hist = matRad_daoVec2ApertureInfo_bixWeightOnly(apertureInfo_hist,apertureInfo_hist.apertureVector);
    
    % get dose for this particular trajectory
    d_hist = matRad_backProjection(apertureInfo_hist.bixelWeights,dij,options);
    
    % update dSum and d2Sum
    dSum    = dSum + d_hist;
    d2Sum   = d2Sum + d_hist.^2;
    
end

% calculate mean and variance of dose
dMean_MC   = dSum./nHistories;
dVar_MC    = (d2Sum-nHistories.*dMean_MC.^2)./(nHistories-1);

end

