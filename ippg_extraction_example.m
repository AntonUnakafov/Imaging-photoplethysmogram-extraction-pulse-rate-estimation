% This script demonstrates the steps of the iPPG extraction and processing by 
% "Imaging photoplethysmogram extraction and pulse rate estimation" package
% developed according to the framework in [a]
%
% CITING THE CODE
% [a] Unakafov AM. Pulse rate estimation using imaging photoplethysmography: generic framework and comparison of methods on a publicly available dataset. Biomedical Physics & Engineering Express. 2018;4(4):045001.
% 
%

% To use this script please either download the DEAP dataset 
% from http://www.eecs.qmul.ac.uk/mmv/datasets/deap/ 
% or use your own video file. 
% In the latter case please change the settings respectively

fontSize = 10;
fontType = 'Arial';
lineWidth = 2;

%% Step 1: extracting raw color signals from the video 

% please see [a] for the description of HSV and STD masking
% here the settings from [a] are used
videoSettings = struct('manualROIselection',  true, ...
                       'isSTDmaskingOn',  true, ...
                       'stdCoef', 1.5, ... 
                       'isHSVmaskingOn',  true, ...
                       'hsvMin', [0.00, 0.09, 0.34], ...
                       'hsvMax', [0.13, 0.52, 1.00]);

%Please set a correct path to the video file: 
video_path = 'C:\Documents\DEAP\Video\s01\s01_trial01.avi';
rawColorSignal = extract_color_channels_from_video(video_path, videoSettings);

%% Steps 2-4: compute and process imaging photoplethysmogramm (iPPG) 
VIDEO_SR = 50.0; % video sampling rate
if (VIDEO_SR > 50)
    FFT_WINDOW_SIZE = 2048; 
elseif (VIDEO_SR > 25)
    FFT_WINDOW_SIZE = 1024; 
else
    FFT_WINDOW_SIZE = 512; 
end    
ippgSettings = struct('PROCESSING', struct('SPA_DETRENDING', 2^0, ...
                                           'PRE_MA_FILTER', 2^1, ...
                                           'PRE_BANDPASS_FILTER', 2^2, ...
                                           'POST_MA_FILTER', 2^3, ...
                                           'POST_BANDPASS_FILTER', 2^4, ...
                                           'OUTLIERS_SUPPRESSION', 2^5, ...
                                           'POST_WAVELET_FILTER', 2^6, ...                                                
                                           'ADAPTIVE_BAND_FILTER', 2^7), ...
                      'EXTRACTION', struct('GREEN', 1, ...
                                           'ICA', 2, ...
                                           'CHROM', 3, ...
                                           'POS', 4, ...
                                           'G_MINUS_R', 5, ...
                                           'AGRD', 6), ...                                                 
                   'minFreq', 0.65, ...
                   'maxFreq', 4.00,...                    
                   'samplingRate', VIDEO_SR, ...                    
                   'processing', 0, ...
                   'extractionMethod', 0, ...  
                   'maFilterLength', 5, ...
                   'bandpassFilterOrder', 127, ... 
                   'stdWindowLength', fix(VIDEO_SR), ...
                   'waveletFilterWidth', [2, 5], ... 
                   'waveletFilterAveragingLength', 20*VIDEO_SR, ...                 
                   'fftWindow', FFT_WINDOW_SIZE, ...
                   'fftShiftSize', ceil(FFT_WINDOW_SIZE/2));

% 1. BASIC PARAMETERS
% 1.1 minFreq - minimal expected pulse rate (0.5-0.7 Hz for humans);
% 1.2 maxFreq - maximal expected pulse rate (4.0 Hz for humans);
% 1.3 samplingRate - frame rate of the video amd of iPPG signal (in Hz)
%
% 2. PROCESSING& IPPG EXTRACTION PARAMETERS
% 2.1 extractionMethod - method for iPPG extraction, please use one of the 
%     constants defined in ippgSettings.EXTRACTION, see below
% 2.2 processing - methods for color signals and iPPG processing, please use
%     a sum of the constants defined in ippgSettings.PROCESSING, see below
% 2.3 maFilterLength - length of moving average filter (if used)
% 2.4 bandpassFilterOrder - order of bandpass filter (if used)
%     for bandpassFilterOrder > 63, (bandpassFilterOrder+1)-point FIR filter is used
%     for bandpassFilterOrder <= 63 IIR Butterworth filter is used
% 2.5 stdWindowLength - length for a window of artefact rejection 
%     (points with amplitude above 3*std are set to 3*std, 
%     std is computed in window of stdWindowLength)
% 2.6 waveletFilterWidth - width of wavelet filter (if used), please
%     see [a] for details
% 2.7 waveletFilterAveragingLength - averaging interval of wavelet filter
%     (if used), please see [a] for details
%
% pulse rate estimation paremeters:
% 3.1. fftWindow: intereger length of the window (in samples) used for 
%      pulse rate estimation via Discrete Fourier Transform (DFT). 
%      Resolution of pulse rate estimation is 60*VIDEO_SR/fftWindow BPM,  
%      thus taking too low fftWindow results in a crude estimation. 
%      Please also avoid too large fftWindow which results in low temporal resoultion 
%      (averaging over too long intervals)
%
% 3.2. fftShiftSize: integer number of samples determining shift of  the sliding 
%      window for computing pulse rate estimates. The larger shift the less 
%      overlapped the windows are. Set fftShiftSize < fftWindow;
%
%     - fftWindow - length of the window  for computing
%       average pulse rate
%     - fftShiftSize - the number of samples determining shift of 
%       the sliding window for computing pulse rate estimates

ippgSettings.processing = ippgSettings.PROCESSING.OUTLIERS_SUPPRESSION + ...                              
                          ippgSettings.PROCESSING.POST_BANDPASS_FILTER;              
ippgSettings.extractionMethod = ippgSettings.EXTRACTION.POS;                      
iPPG = compute_ippg(rawColorSignal, ippgSettings);

xValues = (1:length(iPPG))/ippgSettings.samplingRate;

figure
set( axes,'fontsize', fontSize,  'FontName', fontType);%'FontName', 'Times');
plot(xValues, iPPG, 'LineWidth', lineWidth)
set( gca, 'fontsize', fontSize, 'FontName', fontType);
xlabel( ' time [s] ', 'fontsize', fontSize, 'FontName', fontType);
ylabel( ' imaging photoplethysmogram [a.u.] ', 'fontsize', fontSize, 'FontName', fontType);
title('iPPG extracted using POS method', 'fontsize', fontSize, 'FontName',fontType)
axis tight

%% Steps 5: pulse rate estimation 
% pulse rate estimation using Discrete Fourier Transform:
[hrEstimatedFFT, spectrum] = DFT_pulse_rate_estimate(iPPG, ippgSettings);
% pulse rate estimation using wavelet transform:
[hrEstimatedWavelet, hrMomentaryWavelet] = wavelet_pulse_rate_estimate(iPPG, ippgSettings);

xValues = (1:length(hrMomentaryWavelet))/ippgSettings.samplingRate;
figure
set( axes,'fontsize', fontSize,  'FontName', fontType);%'FontName', 'Times');
plot(xValues, hrMomentaryWavelet, 'LineWidth', 1)
set( gca, 'fontsize', fontSize, 'FontName', fontType);
xlabel( ' time [s] ', 'fontsize', fontSize, 'FontName', fontType);
ylabel( ' pulse rate [BPM] ', 'fontsize', fontSize, 'FontName', fontType);
title('Wavelet-based pulse rate estimation', 'fontsize', fontSize, 'FontName',fontType)
axis tight
