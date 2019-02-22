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

% compute fft of sample signal
[fft_abs_sample,~,fft_freq] = matRad_FFT(data,options.doWindowing);

% set up arrays
fft_absM_simulated      = zeros(numel(fft_freq),1);
fft_abs2M_simulated     = zeros(numel(fft_freq),1);

PijM_simulated      = zeros(nSubPhases,nSubPhases,1);
Pij2M_simulated     = zeros(nSubPhases,nSubPhases,1);

PiM_simulated   = zeros(nSubPhases,1);
Pi2M_simulated  = zeros(nSubPhases,1);

histT_simulatedArray    = zeros(numel(binEdgesMin),nHistories);

% set up simulated data
data_simulated = data;




%%

% simulate a sequence of the Markov chain with Monte Carlo
l_simulated = matRad_runMarkovChain(Pij,nSteps,l_init,false);
data_simulated.l_sample = l_simulated;

% figure out how many long run should be
nStepsPerRun = 2^nextpow2(10./data.deltaT_sample);

% figure out how many runs there are
nRuns = floor(nSteps./nStepsPerRun);

% reproduce data structures
data_cut = data;
data_simulated_cut = data_simulated;

% initialize means
fft_absM_sample_cut     = zeros(nStepsPerRun/2+1,1);
fft_abs2M_sample_cut    = zeros(nStepsPerRun/2+1,1);
fft_absM_simulated_cut  = zeros(nStepsPerRun/2+1,1);
fft_abs2M_simulated_cut = zeros(nStepsPerRun/2+1,1);

% loop over runs
for run = 1:nRuns
    
    % determine run indices
    runIndices = (1:nStepsPerRun) + (run-1).*nStepsPerRun;
    
    % make cuts
    data_cut.l_sample = data.l_sample(runIndices);
    data_simulated_cut.l_sample = data_simulated.l_sample(runIndices);
    
    % do FFT
    [fft_abs_sample_cut,~,fft_freq_cut] = matRad_FFT(data_cut,options.doWindowing);
    [fft_abs_simulated_cut,~,~] = matRad_FFT(data_simulated_cut,options.doWindowing);
    
    % add in means
    fft_absM_sample_cut     = fft_absM_sample_cut+fft_abs_sample_cut./nRuns;
    fft_abs2M_sample_cut    = fft_abs2M_sample_cut+fft_abs_sample_cut.^2./nRuns;
    fft_absM_simulated_cut  = fft_absM_simulated_cut+fft_abs_simulated_cut./nRuns;
    fft_abs2M_simulated_cut = fft_abs2M_simulated_cut+fft_abs_simulated_cut.^2./nRuns;
    
end

% determine std
fft_absS_sample_cut     = (nRuns./(nRuns-1)).*(fft_abs2M_sample_cut-fft_absM_sample_cut.^2);
fft_absS_simulated_cut  = (nRuns./(nRuns-1)).*(fft_abs2M_simulated_cut-fft_absM_simulated_cut.^2);

for history = 1:nHistories
    
    % simulate a sequence of the Markov chain with Monte Carlo
    l_simulated = matRad_runMarkovChain(Pij,nSteps,l_init,false);
    data_simulated.l_sample = l_simulated;
    
    %{
    % extract probability matrices from the simulated sequence
    model_simulated = matRad_generateProbMat(data_simulated);
    
    % compute mean, mean of square
    PijM_simulated  = PijM_simulated+model_simulated.Pij_deltaTSample./nHistories;
    Pij2M_simulated = Pij2M_simulated+model_simulated.Pij_deltaTSample.^2./nHistories;
    
    PiM_simulated   = PiM_simulated+model_simulated.Pi_deltaTSample./nHistories;
    Pi2M_simulated  = Pi2M_simulated+model_simulated.Pi_deltaTSample.^2./nHistories;
    
    
    % compute histogram of periods
    histT_simulatedArray(:,history) = matRad_periodHist(l_simulated,data_simulated.deltaT_sample,data_simulated.indices.nPosSubPhases/2,binEdgesMin,binWidths);
    %}
    
    % compue fft of simulated signal
    [fft_abs_simulated,~,~] = matRad_FFT(data_simulated,options.doWindowing);
    
    % insert in arrays
    fft_absM_simulated  = fft_absM_simulated+fft_abs_simulated./nHistories;
    fft_abs2M_simulated = fft_abs2M_simulated+fft_abs_simulated.^2./nHistories;
    
    matRad_progress(history,nHistories);
end

% calculate mean and standard deviation of simulated matrices
Pij_mean = PijM_simulated;
Pij_std = (nHistories./(nHistories-1)).*(Pij2M_simulated-PijM_simulated.^2);
Pij_mean_std = Pij_std./sqrt(nHistories);
%Pij_mu4 = sum((Pij_simulatedArray-repmat(Pij_mean,1,1,nHistories)).^4,3)./(nHistories-1);
%Pij_std_std = sqrt((Pij_mu4-(nHistories-3).*(Pij_std.^4)./(nHistories-1))./nHistories);

Pi_mean = PiM_simulated;
Pi_std = (nHistories./(nHistories-1)).*(Pi2M_simulated-PiM_simulated);

% insert variables in struct
model.Pij_std   = Pij_std;

end

