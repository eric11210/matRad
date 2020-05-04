function model = matRad_generateMotionModel(options)

switch options.data.origin
    case 'file'
        % read motion data from file
        fprintf('matRad: Reading motion data... ');
        data = matRad_readMotionData(options.data.fileInfo,options.processing.percExtTarg);
        
    case 'function'
        % generate data from function
        fprintf('matRad: Generating motion data... ');
        data = matRad_generateMotionData(options.data);
end
fprintf('Done!\n')

% process data, separating into training and testing data
fprintf('matRad: Processing motion data... ');
[data_train,data_test,parametersNoGood] = matRad_processMotionData(data,options.processing,options.data.fileInfo);
fprintf('Done!\n')

% return if parametersNoGood
if parametersNoGood
    fprintf('matRad: nTimeFracs too large for the choice of Markov frequency!');
    % insert outputs
    if options.hist.doHist
        model.pSum          = NaN;
        model.p             = NaN(size(options.hist.timePoints));
        model.chiSquares    = NaN(size(options.hist.timePoints));
    end
    if options.convTime.doConvTime
        model.convergeT = NaN;
    end
    return
end

% generate matrices using training data
fprintf('matRad: Generating probability matrix... ');
model = matRad_generateProbMat(data_train);
fprintf('Done!\n')

% clean up unused matrix elements
fprintf('matRad: Cleaning up probability matrices and indices... ');
model = matRad_cleanProbMat(model);
fprintf('Done!\n');

if options.convTime.doConvTime
    % determine time to convergence
    fprintf('matRad: Determing convergence time... ');
    model = matRad_motionConvTime(model,options.convTime);
    fprintf('Done!\n')
end

if options.hist.doHist
    % calculate predicted and observed position histograms (as a function of t)
    % this time use testing data
    fprintf('matRad: Calculating position histograms... ');
    model = matRad_calcPosHist(model,data_test,options.hist);
    fprintf('Done!\n')
end

if options.FFT.doFFT
    % calculate the Fourier transform of the testing data
    [model.FFT.fft_abs,model.FFT.fft_phase,model.FFT.fft_freq] = matRad_FFT(data_test,options.FFT.doWindowing);
end

% insert options in model struct
model.options = options;

%{
% estimate standard deviation of matrices
fprintf('matRad: Estimating standard deviation... ');
model = matRad_probMatStd(model,data,options.stdANDfft);
fprintf('Done!\n')
%}
end

