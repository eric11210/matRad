function [fft_abs,fft_phase,fft_freq] = matRad_FFT(signal,samplingFrequency,doDecimation,decFactor)

if nargin < 3
    doDecimation = true;
end

% get the next power of 2 from the original signal length to increase
% efficiency of FFT
n = 2^nextpow2(numel(signal));
%n = numel(signal);

% get Nyquist frequency
fNyq = samplingFrequency./2;

% get frequency resolution
fRes = samplingFrequency./n;

% get the frequencies of each bin
fft_freq = 0:fRes:fNyq;

if nargin < 4
    decFactor = floor(0.05./fRes); 
end

% get the raw FFT
% this is a complex-valued function
fft_raw = fft(signal,n);

% assuming the signal is real, this is symmetrical (technically hermitian)
% therefore double the absolute values of positive frequencies
fft_raw             = fft_raw./n;
fft_raw             = fft_raw(1:(n/2+1));
fft_raw(2:(n/2))    = 2.*fft_raw(2:(n/2));

% get the absolute value
fft_abs             = abs(fft_raw);

% get the phase angle
fft_phase = angle(fft_raw);


% maybe do decimation in time-domain? (window/padding)
if doDecimation
    % do R-decimation described in Ernst's thesis
    
    % first find indices we are interested in
    if mod(decFactor,2)
        decInd1 = ceil(decFactor./2);
    else
        decInd1 = decFactor./2+1;
    end
    decInd = decInd1:decFactor:(n/2+1);
    
    % now do moving average
    fft_abs     = movmean(fft_abs,decFactor);
    fft_phase   = movmean(fft_phase,decFactor);
    
    % now pick out the elements we want to pick, i.e. every decFactor one
    fft_abs     = fft_abs(decInd);
    fft_phase   = fft_phase(decInd);
    fft_freq    = fft_freq(decInd);
    
end

% set to 0 all phases with magnitude below threshold
th = 1e-2;
fft_phase(fft_abs < th) = 0;

end

