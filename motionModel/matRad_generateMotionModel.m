function model = matRad_generateMotionModel(options)

% read motion data from file
fprintf('matRad: Reading motion data... ');
data = matRad_readMotionData(options.fileInfo);
fprintf('Done!\n')

% process data
fprintf('matRad: Processing motion data... ');
data = matRad_processMotionData(data,options.processing);
fprintf('Done!\n')

% generate matrices
fprintf('matRad: Generating probability matrix... ');
model = matRad_generateProbMat(data);
fprintf('Done!\n')

% estimate standard deviation of matrices
fprintf('matRad: Estimating standard deviation... ');
model = matRad_probMatStd(model,data,options.std);
fprintf('Done!\n')

% determine time to convergence
fprintf('matRad: Determing convergence time... ');
model = matRad_motionConvTime(model,options.convTime);
fprintf('Done!\n')

end

