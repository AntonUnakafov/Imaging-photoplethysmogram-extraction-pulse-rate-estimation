% @brief auxilary function wavelet_init_scales generates scales to be used
% for wavelet filtering of iPPG signal or wavelet-based pulse rate estimation
%
% INPUT
%   - ippgSettings - structure with filtering settings and signal properties. 
%     It must contain the following fields:
%       - samplingRate - sampling rate of iPPG signal in Hz,
%       - minFreq - minimal expected pulse rate (0.5-0.7 for humans)
% OUTPUT:
%   - scales - wavelet scales for iPPG signal processing
%
function scales = wavelet_init_scales(ippgSettings)
  MorletFourierFactor = 4*pi/(6+sqrt(2+6^2)); 
  sc0 = 1/(0.5*ippgSettings.samplingRate*MorletFourierFactor); %take default smallest scale (corresponding to Nyquist frequency)
  ds = 0.03125; %32 scales per octave
  scMax = 1/(0.5*ippgSettings.minFreq*MorletFourierFactor); % we do not consider low freq components
  nSc =  fix(log2(scMax/sc0)/ds); 
  scales = {sc0, ds, nSc}; % we use default formula for scales: sc0*2.^((0:nSc-1)*ds)
end