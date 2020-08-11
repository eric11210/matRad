function [pMU_MC,MUMean_MC,MUVar_MC] = matRad_MUMC(apertureInfo,nHistories,p)

% shuffle rng
rng('shuffle');

% default percentiles
if nargin < 6
    p = [5 25 50 75 95];
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

% determine the times for each step, which are the same for all histories
tSimulated_hist = [0; cumsum(tTrans)];

% calculate the number of steps to take (it's possible that there is really
% one more step than necessary, but it's good to have more steps than we
% require than not enough)
nSteps = 1+ceil(tSimulated_hist(end)./apertureInfo.motionModel.deltaT_sample);

% set up histogram of total MU
MU_hist = zeros(nHistories,1);

%% do MC

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
    
    % calculate MU
    for i = 1:numel(apertureInfo_hist.beam)
        MU_hist(hist) = MU_hist(hist)+apertureInfo_hist.beam(i).shape{lSimulated_hist(i)}.MU;
    end
    
end

% calculate mean, variance, and percentiles of MU
MUMean_MC   = sum(MU_hist)./(nHistories);
MUVar_MC    = (sum(MU_hist.^2)-nHistories.*MUMean_MC.^2)./(nHistories-1);
pMU_MC      = prctile(MU_hist,p);

end

