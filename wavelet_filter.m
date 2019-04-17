% @brief wavelet_filter filters the iPPG signal using a two-step wavelet 
% filter suggested in [1]. The filter settings and properties of the signal 
% are passed in structure ippgSettings
%
% INPUT
%   - signal - iPPG signal as a row or column vector;
%   - ippgSettings - structure with filtering settings and signal properties. 
%     It must contain the following fields:
%       - samplingRate - sampling rate of iPPG signal in Hz,
%       - minFreq - minimal expected pulse rate (0.5-0.7 for humans)
%       - waveletFilterWidth - 1x2 array, containing widths of the global (1) 
%         and local (2) wavelet filters, see [1] for details. 
%         Originally, it was suggested to use waveletFilterWidth = [1,4]
%       - waveletFilterAveragingLength - length of averaging used by global
%         wavelet filter. Originally, waveletFilterAveragingLength = 15*samplingRate
%         was suggested
% OUTPUT:
%   - filteredSignal - iPPG undergone wavelet filtering
%
% REFERENCES:
% [1] Bousefsaf F, Maaoui C, Pruski A (2016) 
% Peripheral vasomotor activity assessment using a continuous 
% wavelet analysis on webcam photoplethysmographic signals. 
% Bio-medical materials and engineering. 1;27(5):527-38.
%
function filteredSignal = wavelet_filter(signal, ippgSettings) 
  % prepare parameters for the wavelet transform
  sig = struct('val',signal, 'period', 1/ippgSettings.samplingRate); 
  waveName = {'morl', []};
  sca = wavelet_init_scales(ippgSettings);
  cwtstruct = cwtft(sig, 'wavelet', waveName, 'scales', sca);
  %global filtering
  filteredWaveletCoeff = adaptive_band_gaussian_filter(cwtstruct.cfs, ...
                      ippgSettings.waveletFilterWidth(1), ippgSettings.waveletFilterAveragingLength);
  %local filtering                    
  cwtstruct.cfs  = adaptive_band_gaussian_filter(filteredWaveletCoeff, ...
                      ippgSettings.waveletFilterWidth(2), 1); 
  filteredSignal = icwtft(cwtstruct);  
end  

function waveletCoef = adaptive_band_gaussian_filter(waveletCoef, a, averagingPeriod)
  energyProfile = abs(real(waveletCoef.^2)); 
 
  if (averagingPeriod > 1)
    %filter energyProfile with moving average
    energyProfileMean = movmean(energyProfile, averagingPeriod, 2, 'Endpoints','discard'); %2 stands for average along time
    startPart = energyProfileMean(:, 1);
    energyProfileMean = [startPart(:,ones(1,averagingPeriod-1)), energyProfileMean];
  else
    energyProfileMean = energyProfile;
  end  
  [~, maxScale] = max(energyProfileMean); %maximal scale for each moment of time t
  
  [nScales, nTimeSteps] = size(energyProfileMean);  
  gaussFilter = zeros(nScales, nTimeSteps);  
  
  halfWindowLength = round(maxScale/2);

  % filter wavelet coefficient using normal Gaussian window
  minFilterIndex = 1 - maxScale;
  maxFilterIndex = nScales - maxScale;
  
  for iTime = 1:nTimeSteps   
    x = (minFilterIndex(iTime):maxFilterIndex(iTime))/halfWindowLength(iTime);
    gaussianWindow = exp(-0.5*((a*x).^2));
    gaussFilter(:, iTime) = gaussianWindow;
  end  
  waveletCoef = gaussFilter.*waveletCoef;
end  




