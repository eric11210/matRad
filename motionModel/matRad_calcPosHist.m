function model = matRad_calcPosHist(model,data,options)

%% initialize

% determine integer number of time steps to take
numSteps = round(options.timePoints./model.deltaT_sample);
% determine actual time points
roundTimePoints = numSteps.*model.deltaT_sample;

%% determine initial distribution

% find position phase with max probability
[~,initPhase] = max(accumarray(data.indices.subPhase2PosPhase,model.Pi_deltaTSample));

% initial distribution is all subphases giving the initial phase equally likely
initInd = data.indices.subPhase2PosPhase == initPhase;
initDist = sum(model.Pi_deltaTSample(:,:,1),2);
initDist(~initInd) = 0;
initDist = initDist./sum(initDist);

%% calculated Markov chain probabilities

% initialize probabilities
prob_pred = zeros(data.indices.nPosPhases,numel(numSteps));

% loop through time points
for i = 1:numel(numSteps)
    
    % determine number of steps
    n = numSteps(i);
    
    % calculate distribution across all subphases
    distSubPhase = initDist'*model.Pij_deltaTSample^n;
    
    % sum over subphases
    distPhase = accumarray(data.indices.subPhase2PosPhase,distSubPhase);
    
    % insert into histogram
    prob_pred(:,i) = distPhase;
    
end

%% observed trace: calculate histograms, chi square

% convert l_sample to p_sample
p_sample = data.indices.subPhase2PosPhase(data.l_sample);

% calculate histograms
[hist_obs,initDists] = calcHist(p_sample,data.l_sample,initPhase,data.indices.nPosPhases,data.indices.nSubPhases,numSteps);

% determine last non-zero entry
lastStep = find(sum(hist_obs,1) == 0,1)-1;
if isempty(lastStep)
    lastStep = numel(numSteps);
end

% calculate matrix exponentials
Pijn_deltaTSample = zeros(data.indices.nSubPhases,data.indices.nSubPhases,lastStep);
for i = 1:lastStep
    
    % determine number of steps
    n = numSteps(i);
    
    Pijn_deltaTSample(:,:,i) = (model.Pij_deltaTSample^n);
end

% now calculate predicted histograms
hist_pred       = repmat(sum(hist_obs,1),data.indices.nPosPhases,1).*prob_pred;
hist_var_pred   = repmat(sum(hist_obs,1),data.indices.nPosPhases,1).*prob_pred.*(1-prob_pred);

% now calc chi squares
chiSquares = calcChiSquares(hist_obs,hist_pred,hist_var_pred);


%% Markov chain Monte Carlo trace: calculate histograms, chi square

nHistories = 500;

% initialize chi squares
chiSquaresMC = zeros(nHistories,numel(numSteps));

for history = 1:nHistories
    
    % do MC
    l_simulated = matRad_runMarkovChain(model.Pij_deltaTSample,numel(p_sample),data.l_sample(1),false);
    
    % convert l_sample to p_sample
    p_MCsample = data.indices.subPhase2PosPhase(l_simulated);
    
    % calculate histograms
    [hist_MCobs,initMCDists] = calcHist(p_MCsample,l_simulated,initPhase,data.indices.nPosPhases,data.indices.nSubPhases,numSteps);
    
    
    % for each history, calculate chi square
    
    % start by calculating predicted histograms
    hist_pred       = repmat(sum(hist_MCobs,1),data.indices.nPosPhases,1).*prob_pred;
    hist_var_pred   = repmat(sum(hist_MCobs,1),data.indices.nPosPhases,1).*prob_pred.*(1-prob_pred);
    
    % now calc chi squares
    chiSquaresMC(history,:) = calcChiSquares(hist_MCobs,hist_pred,hist_var_pred);
    
end

%% maximum likelihood Monte Carlo: calculate histograms, chi square

% initialize chi squares
chiSquaresML = zeros(nHistories,numel(numSteps));

% calculate number of observation for each time point
n_obs = sum(hist_obs,1)';

% calculate probabilities for multinomial
pMulti_obs = hist_obs'./repmat(n_obs,1,data.indices.nPosPhases);

for history = 1:nHistories
    
    hist_MLobs = mnrnd(n_obs,pMulti_obs)';
    
    % for each history, calculate chi square
    
    % start by calculating predicted histograms
    hist_pred       = repmat(sum(hist_MLobs,1),data.indices.nPosPhases,1).*prob_pred;
    hist_var_pred   = repmat(sum(hist_MLobs,1),data.indices.nPosPhases,1).*prob_pred.*(1-prob_pred);
    
    % now calc chi squares
    chiSquaresML(history,:) = calcChiSquares(hist_MLobs,hist_pred,hist_var_pred);
    
end


%% calculate p values

% initialize vector
p = zeros(1,numel(numSteps));

% loop through time points
for i = 1:lastStep
    
    % p is the probability of getting a worse disagreement than that
    % observed, just by random chance (under the assumption that our model
    % is correct)j
    % use our MC calculated chi squares instead of the chi square
    % distribution (ours probably doesn't follow chi square)
    p(i) = nnz(chiSquaresMC(:,i) > chiSquares(i))./nHistories;
    
end

pSum = nnz(sum(chiSquaresMC,2) > sum(chiSquares))./nHistories;

%% put hists into structure

model.chiSquares    = chiSquares;
model.chiSquaresMC  = chiSquaresMC;
model.chiSquaresML  = chiSquaresML;
model.p             = p;
model.pSum          = pSum;


end

%% helper functions

% calculate chi square
function chiSquares = calcChiSquares(hist_obs,hist_exp,hist_var_exp)

% calc square of diff over expected
sqDiffByPred = ((hist_obs-hist_exp).^2)./hist_var_exp;

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


