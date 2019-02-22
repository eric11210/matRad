function [fft_abs,fft_phase,fft_freq] = matRad_FFT(data,doWindowing)

% default is to do windowing
if nargin < 3
    doWindowing = true;
end

% extract the signal from data struct, convert to position, not phase
%signal = data.indices.posSubPhase2Pos(data.indices.subPhase2PosSubPhase(data.l_sample));
signal = data.indices.posPhase2Pos(data.indices.subPhase2PosPhase(data.l_sample));

if doWindowing
    % window the signal before FFT to get better results
    % use a better kernel than default sinc
    
    signal = signal.*hamming(numel(signal),'periodic');
end

% get the next power of 2 from the original signal length to increase
% efficiency of FFT
n = 2^nextpow2(numel(signal));

% get the sampling frequency from data struct
samplingFrequency = 1./data.deltaT_sample;

% get Nyquist frequency
fNyq = samplingFrequency./2;

% get frequency resolution
fRes = samplingFrequency./n;

% get the frequencies of each bin
fft_freq = 0:fRes:fNyq;

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

% set to 0 all phases with magnitude below threshold
th = 1e-2;
fft_phase(fft_abs < th) = 0;

end

