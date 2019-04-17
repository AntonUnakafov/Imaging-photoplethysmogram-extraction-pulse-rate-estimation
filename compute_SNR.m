% @brief compute_SNR computes Signal-to-Noise Ration (SNR) of the
% ppg signal given the true pulse rate (only pulse rate components 
% assumed to be the signal of interest). See [1] for details
%
% INPUT
%   - spectrum - cell array 1 x nM of signal spectrums computed by 
%     nM various methods of ppg extraction or spectrum estimation
%     each cell contains a matrix nF x nT representing the spectrum with nF 
%     frequences and nT temporal points     
%   - freq - 1 x nF array of frequences represented in spectrums
%   - trueHR - 1 x nT array of true pulse rate values
%   - nBins - number of bins around the bin corresponding the true pulse
%     rate value, that are counted as signal and not as noise
%
% OUTPUT:
%   - matrix nT x nM, characterizing dynamics of SNR for each of the
%     methods
%
% REFERENCES:
% [1] De Haan G, Jeanne V. (2013)
% Robust pulse rate from chrominance-based rPPG. 
% IEEE Transactions on Biomedical Engineering. 60(10):2878-86.
%
% EXAMPLE OF USE:
%{ 
  % compute trueHR
  % extract channelsIntensity from video
  % define ippgSettings
  settings = ippgSettings;
  for i = 1:ippgSettings.nMethod
    settings.extractionMethod = ippgSettings.extractionMethod(i); 
    settings.processing = ippgSettings.processing(i);
    iPPG = compute_ippg(channelsIntensity, settings);
    [~, spectrum{i}] = DFT_pulse_rate_estimate(iPPG, settings);
  end

  [~, spectrum{i}] = DFT_pulse_rate_estimate(iPPG, ippgSettings);
  [~, spectrum, freq] = computeHRestimate(ppg, ippgSettings)

  %exclude all frequencies outside of heart rate range
  [~, minIndex] = min(abs(freq - ippgSettings.minFreq));
  [~, maxIndex] = min(abs(freq - ippgSettings.maxFreq));
  freq = freq(minIndex:maxIndex);
  spectrum = spectrum(minIndex:maxIndex, :);
  nBins = 4;
  snr = compute_SNR(spectrum, freq, trueHR, nBins)
 %}

function snr = compute_SNR(spectrum, freq, trueHR, nBins)  
  nHRvalue = length(trueHR);
  nMethod = length(spectrum);
  halfNBins = floor(nBins/2);
  
  snr = zeros(nHRvalue, nMethod);
  for iHRvalue = 1:nHRvalue
    [~, hrIndex] = min(abs(freq - trueHR(iHRvalue)));
    [~, secondHarmIndex] = min(abs(freq - 2*trueHR(iHRvalue)));
    
    %take as signal nBins values around the 1st harmonic and 2*nBins+1 around the 2nd:
    signalIndex = [max(hrIndex - halfNBins, 1):(hrIndex + halfNBins), ...
                  (secondHarmIndex - nBins):(secondHarmIndex + nBins)];                       
    signalIndex(signalIndex > length(freq)) = []; %remove too high indices    
    for iMethod = 1:nMethod
      spectrumSquared = spectrum{iMethod}(:, iHRvalue).^2;      
      totalSignalPower = sum(spectrumSquared(signalIndex));
      totalNoisePower = sum(spectrumSquared) - totalSignalPower;
      if (totalNoisePower > 0)
        snr(iHRvalue, iMethod) = 10*log10(totalSignalPower/totalNoisePower);
      else
        snr(iHRvalue, iMethod) = 20;
      end
    end 
  end
end  
 