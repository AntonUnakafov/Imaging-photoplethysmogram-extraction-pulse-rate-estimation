% @brief wavelet_pulse_rate_estimate estimates pulse rate from the iPPG signal 
% using wavelet transform. The function provides two estimates: averaged in
% sliding overlapped windows and momentary pulse rates
%
% INPUT
%   - ippgSignal - iPPG signal as a row vector;
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
%   - instantPulseRate - row vector of non-averaged (momentary) pulse rates
%
function [pulseRate, instantPulseRate] = wavelet_pulse_rate_estimate(ippgSignal, ippgSettings) 
  % wavelet transform
  sig = struct('val',ippgSignal, 'period', 1/ippgSettings.samplingRate); 
  waveName = {'morl', []};
  sca = wavelet_init_scales(ippgSettings);
  cwtstruct = cwtft(sig, 'wavelet', waveName, 'scales', sca);
  
  % determine limits within those we search for pulse rate
  MorletFourierFactor = 4*pi/(6+sqrt(2+6^2));
  freq = 1./(cwtstruct.scales.*MorletFourierFactor);  
  firstScaleIndex = find(freq < ippgSettings.maxFreq, 1, 'first');
  lastScaleIndex = find(freq < ippgSettings.minFreq, 1, 'first'); 

  % compute power spectrum
  energyProfile = abs(real((cwtstruct.cfs.^2)));
  
  % to estimate momentary pulse rate, for every moment of time 
  % we compute frequencies having maximal power within pulse rate limits
  [~, instantPulseRateScales] = max(energyProfile(firstScaleIndex:lastScaleIndex, :));
  instantPulseRate = 60*freq(firstScaleIndex + instantPulseRateScales - 1);
  
  % compute average pulse rates by computing moving average of spectrum
  % and considering only time points with distance ippgSettings.fftShiftSize
  w = ones(1, ippgSettings.fftWindow) ./ ippgSettings.fftWindow;
  energyProfileMean = conv2(energyProfile, w, 'valid');
  [~, maxPos] = size(energyProfileMean);
  pos = 1:ippgSettings.fftShiftSize:maxPos;
  [~, pulseRateScales] = max(energyProfileMean(firstScaleIndex:lastScaleIndex, pos));
  pulseRate = 60*freq(firstScaleIndex + pulseRateScales - 1);  
end  

