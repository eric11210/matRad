p = 6;
f = 19;
m = 1;

data = dlmread(fullfile('ck_processed',sprintf('p%d_f%d_m%d_proc2_pca.csv',p,f,m)));

t = data(:,1);
x = data(:,2);

t_cut = t(57135:end);
x_cut = x(57135:end);

nSubPerPhaseVec = [1:10 20:10:100 200:100:500];
totNumSubPhases = numel(nSubPerPhaseVec);
Pij_meanOfStdVec = zeros(size(nSubPerPhaseVec));
Pij_relMeanOfStdVec = zeros(size(nSubPerPhaseVec));
convergeTVec = zeros(size(nSubPerPhaseVec));

i = 1;
convergeT = 0;
for nSubPerPhase = nSubPerPhaseVec
    
    fprintf('Number of subphases per phase: %d.\n',nSubPerPhase);
    fprintf('Run %d of %d.\n',i,totNumSubPhases);
    
    nPhases = 10;
    
    % convert raw data to phase
    [l_sample,t_sample,deltaT_sample,subPhase2Phase,nSubPhases] = motionData2Phase_NEW(x_cut,t_cut,nPhases,nSubPerPhase);
    
    % construct probability matrix
    [Pij_deltaTSample, Pi_deltaTSample] = generateProbMatrix(l_sample,nSubPhases);
    
    % generate tumour motion model struct
    tmm = generateTmm(Pij_deltaTSample,Pi_deltaTSample,deltaT_sample,subPhase2Phase,nSubPhases);
    
    % estimate standard deviation of probabilities
    nHistories = 10;
    nSteps_sample = numel(l_sample);
    l_init = l_sample(1);
    [Pij_mean,Pij_mean_std,Pij_std,Pij_std_std,Pi_mean,Pi_std] = stdProbMatrix(Pij_deltaTSample,nSubPhases,nHistories,nSteps_sample,l_init);
    Pij_relStd = Pij_std./Pij_mean;
    Pij_relStd(Pij_mean == 0) = 0;
    
    % initial probability vector
    initPVec = zeros(1,nSubPhases);
    initPVec(subPhase2Phase == 1) = 1./nnz(subPhase2Phase == 1);
    
    % determine convergence time to stationary distribution
    targetRelError = 0.01;
    convergeT = convergenceTime(Pij_deltaTSample,Pi_deltaTSample,initPVec,subPhase2Phase,targetRelError,deltaT_sample,convergeT-1); %use last convergence time as initial starting point; convergeT increases with nSubPerPhase
    
    
    % put data in vectors
    Pij_meanOfStdVec(i) = mean(Pij_std(Pij_std ~= 0));
    Pij_relMeanOfStdVec(i) = mean(Pij_relStd(Pij_relStd ~= 0));
    convergeTVec(i) = convergeT;
    i = i+1;
end


save('Convergence times and errors CHANGED.mat')
