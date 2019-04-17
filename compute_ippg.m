% @brief compute_ippg extracts the iPPG signal from the color channels
% by using specified method, pre- and post processing. See [a] for details
%
% CITING THE CODE
% [a] Unakafov AM. Pulse rate estimation using imaging photoplethysmography: generic framework and comparison of methods on a publicly available dataset. Biomedical Physics & Engineering Express. 2018;4(4):045001.
%
% INPUT
%   - rawColorSignal - matrix 3xN containing three color channels (RGB) extracted from video;
%   - ippgSettings - structure with signal properties and iPPG extraction settings. 
%     It must contain the following fields:
%     - PROCESSING: structure defining constants for various procesing settings
%       Must be defined as follows: 
%         PROCESSING = struct('SPA_DETRENDING', 2^0, ...
%                             'PRE_MA_FILTER', 2^1, ...
%                             'PRE_BANDPASS_FILTER', 2^2, ...
%                             'POST_MA_FILTER', 2^3, ...
%                             'POST_BANDPASS_FILTER', 2^4, ...
%                             'OUTLIERS_SUPPRESSION', 2^5, ...
%                             'POST_WAVELET_FILTER', 2^6, ...                                                
%                             'ADAPTIVE_BAND_FILTER', 2^7);
%     - EXTRACTION: structure defining constants for iPPG extraction settings
%       Must be defined as follows:
%         EXTRACTION = struct('GREEN', 1, ...
%                             'ICA', 2, ...
%                             'CHROM', 3, ...
%                             'POS', 4, ...
%                             'G_MINUS_R', 5, ...
%                             'AGRD', 6), ... 
%     - processing - methods for color signals and iPPG processing, please use
%       a sum of the constants defined in ippgSettings.PROCESSING, see below
%     - extractionMethod - method for iPPG extraction, please use one of the 
%       constants defined in ippgSettings.EXTRACTION, see below
%     - samplingRate - sampling rate of iPPG signal in Hz,
%     - minFreq - minimal expected pulse rate (0.5-0.7 Hz for humans)
%     - maxFreq - maximal expected pulse rate (4.0 Hz for humans)
%
% OUTPUT:
%   - finalPPG - vector 1xN of the iPPG values 
%
% REFERENCES: 
% [1] de Haan G, Jeanne V (2013). 
%     Robust pulse rate from chrominance-based rPPG,
%     IEEE Transactions on Biomedical Engineering 60(10): 2878-2886.
% 
% [2] Tarvainen MP, Ranta-Aho PO, Karjalainen PA (2002)
%     An advanced detrending method with application to HRV analysis. 
%     IEEE Transactions on Biomedical Engineering. 49(2):172-5.
%
% [3] McDuff D, Gontarek S, Picard RW (2014). 
%     Improvements in remote cardiopulmonary measurement using a five band digital camera. #
%     IEEE Transactions on Biomedical Engineering 61(10):2593-601. 
%

function finalPPG = compute_ippg(rawColorSignal, ippgSettings)
  dataLength = length(rawColorSignal);
  heartRateBand = [ippgSettings.minFreq ippgSettings.maxFreq]*2/ippgSettings.samplingRate;
  
  %%%%%%%%%%%%%%%%%%%
  %  preprocessing  %
  %%%%%%%%%%%%%%%%%%% 
  if ippgSettings.extractionMethod == ippgSettings.EXTRACTION.AGRD 
    %for adaptive GRD we compute channel weights before preprocessing   
    colorNorm = sqrt(rawColorSignal(1,:).^2 + rawColorSignal(2,:).^2 + rawColorSignal(3,:).^2);
    adaptiveGRDweights = [colorNorm; colorNorm]./(rawColorSignal(1:2,:) + 1);
  end
  
  % detrend signals using Mean-Centering-And-Scaling technique [1]
  % (unconditionally since it's often a required step and it isn't harmful)
  w = ippgSettings.samplingRate;    
  n = conv2(ones(3, dataLength), ones(1, w), 'same');
  meanIntensity = conv2(rawColorSignal, ones(1, w), 'same')./n;
  colorSignal = (rawColorSignal - meanIntensity)./meanIntensity;

  % additional detrending using smoothness prior approach [2]
  if bitand(ippgSettings.processing, ippgSettings.PROCESSING.SPA_DETRENDING)
    % we use crude approximation from the following table to compute lambda
    % SR:    | 30 | 50 | 100
    % lambda | 20 | 50 | 100
    if (ippgSettings.samplingRate < 50) 
      lambda = 1.5*ippgSettings.samplingRate - 25;
    else 
      lambda = ippgSettings.samplingRate;  % for 100 Hz
    end    
    detrendingWindow = 10*ippgSettings.samplingRate;
    colorSignal = smoothness_priors_detrending(colorSignal, lambda, detrendingWindow);
  end
  
  % moving average filter
  if bitand(ippgSettings.processing, ippgSettings.PROCESSING.PRE_MA_FILTER)
    colorSignal = movmean(colorSignal, ippgSettings.maFilterLength, 2);
  end
  
  % band pass filter
  if bitand(ippgSettings.processing, ippgSettings.PROCESSING.PRE_BANDPASS_FILTER) 
    colorSignal = bandpass_filter(colorSignal, heartRateBand, ippgSettings.bandpassFilterOrder);
  end
 
  %%%%%%%%%%%%%%%%%%%%
  %  PPG extraction  %
  %%%%%%%%%%%%%%%%%%%%
  % use single green channel
  if ippgSettings.extractionMethod == ippgSettings.EXTRACTION.GREEN
    iPPG = colorSignal(2, :); 
    
  % use the difference of Green and Red  
  elseif ippgSettings.extractionMethod == ippgSettings.EXTRACTION.G_MINUS_R
    iPPG = colorSignal(2, :) - colorSignal(1, :); 
   
  % use adaptive difference of Green and Red, see [a] for the reference to the original paper
  elseif ippgSettings.extractionMethod == ippgSettings.EXTRACTION.AGRD   
    iPPG = colorSignal(2, :).*adaptiveGRDweights(2, :) - ...
           colorSignal(1, :).*adaptiveGRDweights(1, :); 
  
  % use ICA based extraction [3], see [a] for other references
  elseif ippgSettings.extractionMethod == ippgSettings.EXTRACTION.ICA
    % whitening
    colorSignal = bsxfun(@minus, colorSignal,  mean(colorSignal, 2));
    colorSignal = bsxfun(@rdivide, colorSignal, std(colorSignal, 0, 2) );
    icaCoeff = jadeR(colorSignal);
    icaChannels = icaCoeff*colorSignal;
    iPPG = select_ica_channel(icaChannels, heartRateBand); 
    
  %use CHROM method [1] 
  elseif ippgSettings.extractionMethod == ippgSettings.EXTRACTION.CHROM
    xs = 0.77*colorSignal(1, :) - 0.51*colorSignal(2, :);
    ys = 0.77*colorSignal(1, :) + 0.51*colorSignal(2, :) - 0.77*colorSignal(3, :);         
    alpha = std_ratio_sliding_win(xs, ys, ippgSettings.stdWindowLength);
    iPPG = xs - alpha.*ys;
    
  else %use POS method, see [a] for details and for the reference to the original paper
    xs = colorSignal(2, :) - colorSignal(3, :);
    ys = colorSignal(2, :) + colorSignal(3, :) - 2*colorSignal(1, :);  
    alpha = std_ratio_sliding_win(xs, ys, ippgSettings.stdWindowLength);
    iPPG = xs + alpha.*ys;
  end
  
  %%%%%%%%%%%%%%%%%%%%
  %  postprocessing  %
  %%%%%%%%%%%%%%%%%%%%
  
  refinedPPG = iPPG;
  % moving average filter
  if bitand(ippgSettings.processing, ippgSettings.PROCESSING.POST_MA_FILTER)
    refinedPPG = movmean(refinedPPG, ippgSettings.maFilterLength);
  end

  % band-pass filter
  if (bitand(ippgSettings.processing, ippgSettings.PROCESSING.POST_BANDPASS_FILTER))
    refinedPPG = bandpass_filter(refinedPPG, heartRateBand, ippgSettings.bandpassFilterOrder);
  end
    
  %suppress outliers (unusually high peaks in ppg signal), see [a] for details
  if (bitand(ippgSettings.processing, ippgSettings.PROCESSING.OUTLIERS_SUPPRESSION))
    % make ppg signal zero-mean
    averagingWindow = 2*ippgSettings.samplingRate;
    refinedPPG = refinedPPG - movmean(refinedPPG, averagingWindow);
    %... and with std = 1 
    stdValues = std_sliding_win(refinedPPG, averagingWindow);
    indices = find(stdValues < abs(refinedPPG)/3); 
    stdValues(indices) = abs(refinedPPG(indices))/3;
    refinedPPG = refinedPPG./stdValues;    
  end  
  
  % final processing
  if bitand(ippgSettings.processing, ippgSettings.PROCESSING.POST_WAVELET_FILTER)
    finalPPG = wavelet_filter(refinedPPG, ippgSettings);
  else
    finalPPG = refinedPPG;
  end
  
  %correct sign of final ppg
  finalPPG = finalPPG*calc_ppg_sign(finalPPG);  
end

% select ICA channel corresponding to ppg
function ppgChannel = select_ica_channel(icaChannels, freqBand)
  dataLength = length(icaChannels);
  w = hamming(dataLength)';
  icaChannelWindowed = bsxfun(@times, icaChannels, w);

  %fast Fourier transform (FFT) to obtain the power spectrum
  channelFFT = fft(icaChannelWindowed, [], 2);
  channelSpectrum2side = abs(channelFFT/dataLength);
  channelSpectrum = channelSpectrum2side(:, 1:floor(dataLength/2)+1);
  channelSpectrum(:, 1) = channelSpectrum(:, 1)/2; %divide the 1st coefficiend instead of multiplying other  
  channelSpectrum = bsxfun(@rdivide, channelSpectrum, sum(channelSpectrum, 2));
     
  %find the channel, corresponding to ppg
  freqIndex = floor(1 + freqBand*length(channelSpectrum));  
  maxSpectrumValue = max(channelSpectrum(:, freqIndex(1):freqIndex(2)), [], 2);
  [~, ppgComponentIndex] = max(maxSpectrumValue);
  ppgChannel = icaChannels(ppgComponentIndex, :); 
end

% computes bandpass filter without edge effect
function signalFiltered = bandpass_filter (signal, freqBand, filterOrder)
  if (filterOrder > 11) %use FIR filter    
    b = fir1(filterOrder, freqBand);  
    a = 1;
  else  %use Butterworth filter
    [b, a] = butter(filterOrder, freqBand);
  end
  [nChannel, dataLength] = size(signal);
  signalFiltered = zeros(nChannel, dataLength);
  for iChannel = 1:nChannel
    signalFiltered(iChannel, :) = filtfilt(b, a, signal(iChannel, :));
  end
end

% computes ratio of running standard deviations for two signals
function alpha = std_ratio_sliding_win(x, y, w)
  % n - number of elements in each window
  n = conv(ones(1, length(x)), ones(1, w), 'same');
  
  % compute running variation of x time series
  s = conv(x, ones(1, w), 'same');
  q = x.^ 2;    
  q = conv(q, ones(1, w), 'same');
  varXs = (q - s.^2 ./ n) ./ (n - 1);
    
  % compute running std of y time series  
  s = conv(y, ones(1, w), 'same');
  q = y.^ 2;    
  q = conv(q, ones(1, w), 'same');    
  varYs = (q - s.^2 ./ n)./(n - 1);  
    
  alpha = (varXs./varYs).^ 0.5; %compute std ratio
end

% recovers the correct sign of iPPG signal
% using the fact that in PPG throughs have higher amplitudes than peaks
function ppgSign = calc_ppg_sign(ppg)
  % make the signal zero-mean and normalize its amplitude
  minValue = min(ppg);
  maxValue = max(ppg); 
  meanValue = mean(ppg);   
  ppg = (ppg - meanValue)/(maxValue - minValue);

  ppg(ppg == 0) = 0.00001; % get rid of exact zeros
  signChange = diff(sign(ppg));
  indexUp = find(signChange > 0);
  indexDown = find(signChange < 0);
  if (indexUp(1) > indexDown(1))  % for consistency we start with positive segment
    indexUp = [1, indexUp]; % thus we add first segment if necessary
  end  
  if (indexUp(end) < indexDown(end))  % we also end with positive segment
    indexDown = indexDown(1:end-1); % thus we skip last segment if necessary
  end    
  nSegment = min(length(indexUp), length(indexDown));
  sumMax = 0;
  sumMin = 0;
  for iSegment = 1:nSegment
    sumMax = sumMax + max(ppg(indexUp(iSegment):indexDown(iSegment)) );
    sumMin = sumMin - min(ppg(indexDown(iSegment):indexUp(iSegment+1))) ;    
  end
  if (sumMax > sumMin) % if peaks have higher amplitudes than throughs
    ppgSign = -1;      % signal shold be inverted  
  else
    ppgSign = 1;
  end  
end

