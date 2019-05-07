function [dMean,dVar] = matRad_doseVarianceMC(apertureInfo,dij,nHistories)
% Monte Carlo estimate of the mean and variance of the dose

% initialize sum of dose and square of dose
dSum    = zeros(dij.numOfVoxels,1);
d2Sum   = zeros(dij.numOfVoxels,1);

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
apertureInfo.motionModelOriginal = apertureInfo.motionModel;
apertureInfo = rmfield(apertureInfo,'motionModel');
apertureInfo.motionModel.type = 'single';
apertureInfo.motionModel.numPhases = apertureInfo.numPhases;

% loop over number of histories
for hist = 1:nHistories
    
    % MC simulation of tumour trajectory
    [lSimulated,tSimulated]= matRad_runMarkovChain_Q(apertureInfo.motionModelOriginal.qij,apertureInfo.motionModelOriginal.initProb,sum([apertureInfo.beam.time]));
    
    % insert trajectory in apertureInfo struct
    apertureInfo.motionModel.lSimulated = lSimulated;
    apertureInfo.motionModel.tSimulated = tSimulated;
    
    % get weights for this particular history
    apertureInfo_hist = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfo.apertureVector);
    
    % get dose for this particular trajectory
    d_hist = matRad_backProjection(apertureInfo_hist.bixelWeights,dij,options);
    
    % update dSum and d2Sum
    dSum    = dSum + d_hist;
    d2Sum   = d2Sum + d_hist.^2;
    
end

% calculate mean and variance of dose
dMean   = dSum./nHistories;
dVar    = (d2Sum-nHistories.*dMean.^2)./(nHistories-1);

end
