% @brief smoothness_priors_detrending implements the signal detrending 
% using SMOOTHNESS PRIOR APPROACH (SPA) [1]. 
% Contrary to the original approach, which is defined for short signals,
% here we apply SPA in running windows of given length. The signal is
% detrended in overlapping windows (50% overlap) and then glued.
% To glue overlapping signal parts, hamming window is used.
%
% INPUT
%   - signal - the multi-channel signal for detrending;
%     represented by nC x nT matrix, where nC - number of channels, nT -
%     length of the signal
%   - lambda - control parameter of detrending, see [1] for details;
%   - windowSize - the length of the window for detrending. Depends on the 
%     expected frequences of the trend. Please do not set windowSize to
%     high to avoid computational overload
% OUTPUT:
%   - detrendSignal - detrended signal represented by nC x nT matrix
%
% REFERENCES:
% [1] Tarvainen MP, Ranta-Aho PO, Karjalainen PA (2002)
% An advanced detrending method with application to HRV analysis. 
% IEEE Transactions on Biomedical Engineering. 49(2):172-5.
%

function detrendSignal = smoothness_priors_detrending (signal, lambda, windowSize)
  [nChannels, nSamples] = size(signal);
  % padding of the remaining part of the signal with zeros
  extendedLength = ceil(nSamples/(windowSize/2))*(windowSize/2);
  signal(:, nSamples+1:extendedLength) = zeros(nChannels, extendedLength - nSamples);
    
  I = speye(windowSize);
  D2 = spdiags(ones(windowSize - 2, 1)*[1 -2 1], 0:2, windowSize - 2, windowSize);
  detrendMatrix = (I - inv(I + lambda^2*(D2'*D2)))';
    
  detrendSignal = zeros(nChannels, extendedLength);
  hammingWindow = repmat(hamming(windowSize)', nChannels, 1);
  for j = 1:windowSize/2:extendedLength-windowSize + 1
    % determine the points we are computing now
    points = j:j+windowSize-1;
    % detrend the current signal window
    windowedSignal = signal(:, points)*detrendMatrix;
    % add it to the composite signal
    detrendSignal(:, points) = detrendSignal(:, points) + windowedSignal.*hammingWindow;    
  end  
  
  % reject the padded part of the signal 
  detrendSignal = detrendSignal(:, 1:nSamples);
end
