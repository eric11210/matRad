function model = matRad_calcPosHist(model,data,options)

%% initialize

% determine integer number of time steps to take
numStepsModel_timePoints = round(options.timePoints./model.deltaT_sample);
numStepsData_timePoints = round(options.timePoints./data.deltaT_sample);
numStepsModel_simulate   = round(data.deltaT_sample*numel(data.l_sample)./model.deltaT_sample);
% determine actual time points
roundTimePoints = numStepsModel_timePoints.*model.deltaT_sample;

%% determine initial distribution for histograms

% find position phase with max probability; this is the triggering phase
[~,triggerPhase] = max(accumarray(model.indices.subPhase2PosPhase,model.Pi_deltaTSample));

% initial distribution is all subphases giving the initial phase equally likely
triggerInd = model.indices.subPhase2PosPhase == triggerPhase;
triggerDist = model.Pi_deltaTSample;
triggerDist(~triggerInd) = 0;
triggerDist = triggerDist./sum(triggerDist);

% determine distribution of initial subphases for simulation
initPhase = data.indices.subPhase2PosPhase(data.l_sample(1));
initInd = model.indices.subPhase2PosPhase == initPhase;
initDist = model.Pi_deltaTSample;
initDist(~initInd) = 0;
initDist = initDist./sum(initDist);

%% Markov chain Monte Carlo trace: expected number of visits to triggering phase

nHistories = 1000;

% initialize number of triggering phase
numTriggerPhaseMC  = zeros(1,numel(numStepsModel_timePoints));

% now calculate cumulative distribution
cumInitSimDist = cumsum(initDist);

for history = 1:nHistories
    
    % determine initial subphase for simulation by sampling for the
    % distribution
    r = rand;
    initSubPhase = find(r < cumInitSimDist,1,'first');
    
    % do MC
    l_simulated = matRad_runMarkovChain_P(model.Pij_deltaTSample,numStepsModel_simulate,initSubPhase,false);
    
    % convert l_sample to p_sample
    p_MCsample = model.indices.subPhase2PosPhase(l_simulated);
    
    % loop through time points
    for i = 1:numel(numStepsModel_timePoints)
        
        % determine number of steps
        n = numStepsModel_timePoints(i);
        
        % determine number of times triggering phase occurs
        numTriggerPhaseMC(i) = numTriggerPhaseMC(i)+nnz(p_MCsample(1:(end-n)) == triggerPhase)./nHistories;
    end
end

%% calculate Markov chain histograms

% initialize probabilities
hist_pred = zeros(data.indices.nPosPhases,numel(numStepsModel_timePoints));

% calculate sum of n-step transition matrices for total signal
%polyTransTot        = ones(numel(data.l_sample),1);
%sumNStepTransTot    = polyvalm(polyTransTot,model.Pij_deltaTSample);

% loop through time points
for i = 1:numel(numStepsModel_timePoints)
    
    % determine number of steps
    n = numStepsModel_timePoints(i);
    
    % calculate sum of n-step transition matrices for particular number of
    % steps
    %polyTransPart               = polyTransTot;
    %polyTransPart(1:(end-n))    = 0;
    %sumNStepTransPart           = sumNStepTransTot-polyvalm(polyTransPart,model.Pij_deltaTSample);
    
    % calculate n-step transition matrix for particular number of steps
    nStepTransPart = model.Pij_deltaTSample^n;
    
    % calculate expected number of times to observe the triggering phase
    %histTriggerSubPhase = initDist'*sumNStepTransPart;
    %numTriggerPhase     = sum(histTriggerSubPhase(data.indices.subPhase2PosPhase == triggerPhase));
    
    % calculate probability of observing a phase nSteps after the
    % triggering phase
    distObsSubPhase     = triggerDist'*nStepTransPart;
    distObsPhase        = accumarray(model.indices.subPhase2PosPhase,distObsSubPhase);
    
    % calculate combined histogram: probability of observing a phase nSteps
    % after the triggering phase multiplied by the number of times the
    % triggering phase is observed
    %histPhase = numTriggerPhase.*distObsPhase;
    histPhase = numTriggerPhaseMC(i).*distObsPhase;
    
    % insert into histogram
    hist_pred(:,i) = histPhase;
    
end

%% observed trace: calculate histograms, chi square

% convert l_sample to p_sample
p_sample = data.indices.subPhase2PosPhase(data.l_sample);

% calculate histograms
[hist_obs,initDists] = calcHist(p_sample,data.l_sample,triggerPhase,data.indices.nPosPhases,data.indices.nSubPhases,numStepsData_timePoints);

% scale the histogram of observed (training) data to match the testing data
% this is required due to the difference sampling intervals between the
% training and testing data (former is decimated, latter is not)
hist_obs = hist_obs.*numStepsModel_simulate./numel(data.l_sample);

% now calc chi squares
chiSquares = calcChiSquares(hist_obs,hist_pred);


%% Markov chain Monte Carlo trace: calculate histograms, chi square

nHistories = 1000;

% initialize chi squares
chiSquaresMC        = zeros(nHistories,numel(numStepsModel_timePoints));

for history = 1:nHistories
    
    % determine initial subphase for simulation by sampling for the
    % distribution
    r = rand;
    initSubPhase = find(r < cumInitSimDist,1,'first');
    
    % do MC
    l_simulated = matRad_runMarkovChain_P(model.Pij_deltaTSample,numStepsModel_simulate,initSubPhase,false);
    
    % convert l_sample to p_sample
    p_MCsample = model.indices.subPhase2PosPhase(l_simulated);
    
    % calculate histograms
    [hist_MCobs,initMCDists] = calcHist(p_MCsample,l_simulated,triggerPhase,data.indices.nPosPhases,data.indices.nSubPhases,numStepsModel_timePoints);
    
    % for each history, calculate chi square
    chiSquaresMC(history,:) = calcChiSquares(hist_MCobs,hist_pred);
    
end

%% calculate p values

% initialize vector
p = zeros(1,numel(numStepsModel_timePoints));

% loop through time points
for i = 1:numel(numStepsModel_timePoints)
    
    % p is the probability of getting a worse disagreement than that
    % observed, just by random chance (under the assumption that our model
    % is correct)
    % use our MC calculated chi squares instead of the chi square
    % distribution (ours probably doesn't follow chi square)
    p(i) = nnz(chiSquaresMC(:,i) > chiSquares(i))./nHistories;
    
end

pSum = nnz(sum(chiSquaresMC,2) > sum(chiSquares))./nHistories;

%% put hists into structure

model.chiSquares    = chiSquares;
model.chiSquaresMC  = chiSquaresMC;
model.p             = p;
model.pSum          = pSum;


end

%% helper functions

% calculate chi square
function chiSquares = calcChiSquares(hist_obs,hist_exp)

% calc square of diff over expected
sqDiffByPred = ((hist_obs-hist_exp).^2)./hist_exp;

% set to 0 when expected is zero
sqDiffByPred(hist_exp == 0) = 0;

% sum over bins
chiSquares = sum(sqDiffByPred,1);

end

% calculate G
function Gs = calcGs(hist_obs,hist_exp)

% calc obs times logarithms of ratios
obsMultLnRatios = hist_obs.*ln(hist_obs./hist_exp);

% set to 0 when expected or observed is 0
obsMultLnRatios(hist_obs == 0) = 0;
obsMultLnRatios(hist_exp == 0) = 0;

% sum over bins
Gs = 2.*sum(obsMultLnRatios,1);

end

% calculate observed histograms
function [hist_obs,initDists] = calcHist(b_sample,sb_sample,initBin,nBins,nSubBins,numSteps)

% initialize histograms
hist_obs = zeros(nBins,numel(numSteps));

% initialize initial distribution as func of t
initDists = zeros(nSubBins,numel(numSteps));

% loop through trace
for i = 1:numel(b_sample)
    
    % determine if phase is the initial
    if b_sample(i) == initBin
        
        % initialize steps and indices for this round
        numSteps_i      = numSteps;
        numSteps_i_ind  = 1:numel(numSteps);
        
        % cut steps off if they exceed number of points in the trace
        numSteps_i(numSteps+i > numel(b_sample)) = [];
        numSteps_i_ind(numSteps+i > numel(b_sample)) = [];
        
        % determine phase at numSteps+i
        p_stepped = b_sample(numSteps_i+i);
        
        % determine linear indices
        hist_obs_ind = sub2ind([nBins,numel(numSteps)],p_stepped,numSteps_i_ind');
        
        % insert phases into histogram
        hist_obs(hist_obs_ind) = hist_obs(hist_obs_ind)+1;
        
        % put subphase into histogram
        initDists(sb_sample(i),numSteps_i_ind) = initDists(sb_sample(i),numSteps_i_ind)+1;
        
    end
    
end

end


