options.fileInfo.p          = 6;
options.fileInfo.f          = 19;
options.fileInfo.m          = 1;
options.fileInfo.processed  = true;

options.processing.percExtTarg  = 1;
options.processing.nPosPhases   = 10;
options.processing.nVelPhases   = 10;
options.processing.velBinning   = true;
options.convTime.percRMSD_targ  = 1;
options.std.nHistories          = 10;

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
    
    options.processing.nSubPerVelPhase = 1;
    options.processing.nSubPerPosPhase = 10;
    
    model = matRad_generateMotionModel(options);
    
    % put data in vectors
    Pij_meanOfStdVec(i) = mean(Pij_std(Pij_std ~= 0));
    Pij_relMeanOfStdVec(i) = mean(Pij_relStd(Pij_relStd ~= 0));
    convergeTVec(i) = convergeT;
    i = i+1;
end


save('Convergence times and errors.mat')
