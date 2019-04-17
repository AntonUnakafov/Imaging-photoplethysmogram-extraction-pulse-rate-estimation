% @brief DFT_pulse_rate_estimate estimates pulse rate from the iPPG signal 
% using discrete Fourier transform (DFT). The function provides an estimate
% averaged in sliding overlapped windows
%
% INPUT
%   - ippgSignal - iPPG signal as a row or column vector;
%   - ippgSettings - structure with pulse rate estimation settings and signal properties. 
%     It must contain the following fields:
%       - samplingRate - sampling rate of iPPG signal in Hz,
%       - fftWindow - length of the findow (in samples) for computing
%         average pulse rate
%       - fftShiftSize - the number of samples determining shift of 
%         the sliding window for computing pulse rate estimates
%       - minFreq - minimal expected pulse rate (0.5-0.7 Hz for humans)
%       - maxFreq - maximal expected pulse rate (4.0 Hz for humans)
% OUTPUT:
%   - pulseRate - row vector of pulse rate estimates averaged in sliding (overlapped) windows
%   - spectrum - Matrix of the amplitude spectrum of the signal between 
%     minFreq and maxFreq. The number of columns is the same as in pulseRate,
%     the number of rows is the number of frequencies between minFreq and maxFreq 
%
function [pulseRate, spectrum] = DFT_pulse_rate_estimate(ippgSignal, ippgSettings)
  minHeartFreq = floor(1 + ippgSettings.minFreq*ippgSettings.fftWindow/ippgSettings.samplingRate);
  maxHeartFreq = floor(1 + ippgSettings.maxFreq*ippgSettings.fftWindow/ippgSettings.samplingRate);

  w = hamming(ippgSettings.fftWindow)';  
  dataLength = length(ippgSignal);
  startPos = 1:ippgSettings.fftShiftSize:(dataLength-ippgSettings.fftWindow+1);
  nWin = length(startPos);
  heartFreqIndex = zeros(1, nWin);
  spectrum = zeros(maxHeartFreq - minHeartFreq + 1, nWin);
  
  for iWin = 1:nWin       
    ppgWindowed = ippgSignal(startPos(iWin):startPos(iWin)+ippgSettings.fftWindow-1).*w;   
    ppgFFT = fft(ppgWindowed, [], 2);
    ppgSpectrum2side = abs(ppgFFT/ippgSettings.fftWindow);
    spectrum(:, iWin) = ppgSpectrum2side(minHeartFreq:maxHeartFreq); %take only values in heart rate bandwidth
    %find the maximal frequency (presumably, corresponding to ppg)
    [~, heartFreqIndex(iWin)] = max(spectrum(:, iWin));
  end
  spectrum = bsxfun(@rdivide, spectrum, sum(spectrum));
  % heartFreqIndex + minHeartFreq - 1 is the index of maximal frequency
  % since index = 1 corresponds to frequency 0 Hz (DC), we need to deduce 1: 
  pulseRate = (60*ippgSettings.samplingRate/ippgSettings.fftWindow)*(heartFreqIndex + minHeartFreq - 2);
end