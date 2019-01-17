function model = matRad_probMatStd(model,data,options)

% extract variables from struct
nSubPhases  = model.indices.nSubPhases;
nHistories  = 100;%options.nHistories;
Pij         = model.Pij_deltaTSample;
nSteps      = numel(data.l_sample);
l_init      = data.l_sample(1);

% set up period histograms
binEdgesMin = (0:0.5:10)';
binWidths = 1;
histT_sample = matRad_periodHist(data.l_sample,data.deltaT_sample,data.indices.nPosSubPhases/2,binEdgesMin,binWidths);

% set up arrays
histT_simulatedArray    = zeros(numel(binEdgesMin),nHistories);
Pij_simulatedArray      = zeros(nSubPhases,nSubPhases,nHistories);
Pi_simulatedArray       = zeros(nSubPhases,nHistories);

% set up simulated data
data_simulated = data;

for history = 1:nHistories
    
    % simulate a sequence of the Markov chain with Monte Carlo
    l_simulated = matRad_runMarkovChain(Pij,nSteps,l_init,false);
    
    % extra probability matrices from the simulated sequence
    data_simulated.l_sample = l_simulated;
    model_simulated = matRad_generateProbMat(data_simulated);
    
    % compute histogram of periods
    histT_simulatedArray(:,history) = matRad_periodHist(l_simulated,data_simulated.deltaT_sample,data_simulated.indices.nPosSubPhases/2,binEdgesMin,binWidths);
    
    % insert in arrays
    Pij_simulatedArray(:,:,history) = model_simulated.Pij_deltaTSample;
    Pi_simulatedArray(:,history)    = model_simulated.Pi_deltaTSample;
    
    matRad_progress(history,nHistories);
end

% calculate mean and standard deviation of simulated matrices
Pij_mean = mean(Pij_simulatedArray,3);
Pij_std = std(Pij_simulatedArray,0,3);
Pij_mean_std = Pij_std./sqrt(nHistories);
Pij_mu4 = sum((Pij_simulatedArray-repmat(Pij_mean,1,1,nHistories)).^4,3)./(nHistories-1);
Pij_std_std = sqrt((Pij_mu4-(nHistories-3).*(Pij_std.^4)./(nHistories-1))./nHistories);

Pi_mean = mean(Pi_simulatedArray,2);
Pi_std = std(Pi_simulatedArray,0,2);

% insert variables in struct
model.Pij_std   = Pij_std;

end

