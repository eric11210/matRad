function model = matRad_generateMotionModel(options)

switch options.data.origin
    case 'file'
        % read motion data from file
        fprintf('matRad: Reading motion data... ');
        data = matRad_readMotionData(options.data.fileInfo);
        
    case 'function'
        % generate data from function
        fprintf('matRad: Generating motion data... ');
        data = matRad_generateMotionData(options.data);
end
fprintf('Done!\n')

% process data, separating into training and testing data
fprintf('matRad: Processing motion data... ');
[data_train,data_test] = matRad_processMotionData(data,options.processing);
fprintf('Done!\n')

% generate matrices using training data
fprintf('matRad: Generating probability matrix... ');
model = matRad_generateProbMat(data_train);
fprintf('Done!\n')

% calculate predicted and observed position histograms (as a function of t)
% this time use testing data
fprintf('matRad: Calculating position histograms... ');
model = matRad_calcPosHist(model,data_test,options.hist);
fprintf('Done!\n')

% clean up unused matrix elements
fprintf('matRad: Cleaning up probability matrices and indices... ');
model = matRad_cleanProbMat(model);
fprintf('Done!\n');

% insert options in model struct
model.options = options;

%{
% estimate standard deviation of matrices
fprintf('matRad: Estimating standard deviation... ');
model = matRad_probMatStd(model,data,options.stdANDfft);
fprintf('Done!\n')

% determine time to convergence
fprintf('matRad: Determing convergence time... ');
model = matRad_motionConvTime(model,options.convTime);
fprintf('Done!\n')
%}
end

