% This script demonstrate the performance of the 
% "Imaging photoplethysmogram extraction and pulse rate estimation" package
% developed according to the framework in [a], for the data from [b]
%
% CITING THE CODE
% [a] Unakafov AM. Pulse rate estimation using imaging photoplethysmography: generic framework and comparison of methods on a publicly available dataset. Biomedical Physics & Engineering Express. 2018;4(4):045001.
% [b] Unakafov AM, Moeller S, Kagan I, Gail A, Treue S, Wolf F. Using imaging photoplethysmography for heart rate estimation in non-human primates. PLoS ONE 2018;13(8): e0202581. https://doi.org/10.1371/journal.pone.0202581

%% Initialization
clear all;

global FontSize LineWidth;
FontSize = 9;
LineWidth = 1;

filename = {'Session1', 'Session2', 'Session3', 'Session4', 'Session5', ...
            'Session6', 'Session7', 'Session8', 'Session2IR', 'Session9IR', 'Session10IR'};
nFile = length(filename);

% determines, whether color signals are extracted from an RGB (1) or monochrome IR video (0)
RGB_DATA_FLAG = [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0];
               
% index of the subject (1 - Sun, 2 - Fla, 3 - Mag, 4 - Lin)               
SUBJECT_INDEX = [1, 1, 2, 2, 3, 3, 3, 4, 1, 3, 3];  
               
% index of the session among RGB-sessions of a subject       
SECTION_INDEX_PER_SUBJECT = [1, 2, 1, 2, 1, 2, 3, 1, 0, 0, 0];
   
% frame rate of the video = sampling rate of color signals  
VIDEO_FRAME_RATE = [30, 30, 50, 100, 50, 50, 50, 50, 30, 100, 50];                 
  
% read reference pulse rate values 
REFERENCE_PULSE_RATE = cell(1, nFile);                
for i = 1:nFile
  REFERENCE_PULSE_RATE{i} = csvread(['dataset' filesep filename{i} filesep filename{i} '_REF.csv']);
end


VIDEO_SR = 50.0;     % just to initialize the structure, actual values provided in the loop 
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
                   'fftWindow', 1024, ... 
                   'minFreq', 1.5, ...
                   'maxFreq', 5.0,...                    
                   'samplingRate', VIDEO_SR, ...                    
                   'maFilterLength', 5, ...
                   'bandpassFilterOrder', 127, ... 
                   'stdWindowLength', fix(VIDEO_SR), ...
                   'waveletFilterWidth', [2, 5], ... 
                   'waveletFilterAveragingLength', 20*VIDEO_SR, ...
                   'extractionMethod', 0, ...
                   'processing', 0);
                                   
%number of frames nearest to the fftWindow/2 and corresponding to integer number of seconds in video
DFT_WINDOW_SHIFT = [ 510, 510, 500, 1000, 500, 500, 500, 500, 510, 1000, 500];                     

finalPPG = cell(nFile, 1);

shareErrorBelow3p5BPM = cell(nFile, 1);
shareErrorBelow7BPM = cell(nFile, 1);
corrCoef = cell(nFile, 1);
meanError = cell(nFile, 1);
rmse = cell(nFile, 1);
stdError = cell(nFile, 1);
snr = cell(nFile, 1);
corrPvalue = cell(nFile, 1);

%variables for motion estimation
nBins = {6,8,6,4,6,4,8,6,1,1,1};  % optimal number of bins for computing SNR (selected based on estimation errors)
motionData = cell(nFile, 1);
errorForMotion = cell(nFile, 1);
startPosForMotion = cell(nFile, 1);
endPosForMotion = cell(nFile, 1);
dFreqMotion = cell(nFile, 1);

nSubject = length(unique(SUBJECT_INDEX));
subjectHRtrue = cell(nSubject, 1);
subjectHRestimate = cell(nSubject, 1);
subjectSessionIndex = cell(nSubject, 1);
  
hrTrue = cell(1, nFile);
hrEstimated = cell(1, nFile);
xt = cell(1, nFile);

% estimate pulse rates for the dataset and evaluation of estimates' performance 
for iFile = 1:nFile
  % set iPPG parameters for each file
  if (VIDEO_FRAME_RATE(iFile) <= 50)  
    ippgSettings.fftWindow = 1024;
  else
    ippgSettings.fftWindow = 2048;
  end
  ippgSettings.fftShiftSize = DFT_WINDOW_SHIFT(iFile);
  if (RGB_DATA_FLAG(iFile))  
    ippgSettings.nMethod = 4;
    noWaveletIndex = 1;
    mainMethodIndex = 4;

    ippgSettings.processing = zeros(1, ippgSettings.nMethod, 'uint16');
    ippgSettings.extractionMethod = zeros(1, ippgSettings.nMethod, 'uint16');
  
    ippgSettings.extractionMethod(1) = ippgSettings.EXTRACTION.G_MINUS_R;
    ippgSettings.processing(1) = ippgSettings.PROCESSING.POST_BANDPASS_FILTER + ...
                                 ippgSettings.PROCESSING.OUTLIERS_SUPPRESSION;   
    ippgSettings.extractionMethod(2) = ippgSettings.EXTRACTION.CHROM;
    ippgSettings.processing(2) = ippgSettings.PROCESSING.POST_WAVELET_FILTER + ...
                                 ippgSettings.PROCESSING.OUTLIERS_SUPPRESSION + ...                              
                                 ippgSettings.PROCESSING.POST_BANDPASS_FILTER;
    ippgSettings.extractionMethod(3) = ippgSettings.EXTRACTION.POS;
    ippgSettings.processing(3) = ippgSettings.PROCESSING.POST_WAVELET_FILTER + ...
                                 ippgSettings.PROCESSING.POST_BANDPASS_FILTER + ...
                                 ippgSettings.PROCESSING.OUTLIERS_SUPPRESSION;
    ippgSettings.extractionMethod(4) = ippgSettings.EXTRACTION.G_MINUS_R;
    ippgSettings.processing(4) = ippgSettings.PROCESSING.POST_WAVELET_FILTER + ...
                                 ippgSettings.PROCESSING.POST_BANDPASS_FILTER + ...
                                 ippgSettings.PROCESSING.OUTLIERS_SUPPRESSION;
  else
    ippgSettings.nMethod = 1;   
    noWaveletIndex = 1;
    mainMethodIndex = 1;
    ippgSettings.extractionMethod(1) = ippgSettings.EXTRACTION.GREEN;
    ippgSettings.processing(1) = ippgSettings.PROCESSING.POST_BANDPASS_FILTER + ...
                                 ippgSettings.PROCESSING.OUTLIERS_SUPPRESSION + ...
                                 ippgSettings.PROCESSING.POST_MA_FILTER;
  end  
  ippgSettings.samplingRate = VIDEO_FRAME_RATE(iFile);
  ippgSettings.maFilterLength = ceil(ippgSettings.samplingRate/10);
  
  % read color signals from file
  bFile = fopen(['dataset' filesep filename{iFile} filesep filename{iFile} '.bin']);
  bData = fread(bFile, 'single');
  fclose(bFile);
  nFrames = numel(bData)/3;
  channelsIntensity = reshape(bData, 3, nFrames);
  nAverageHRValues = 1 + floor((nFrames - ippgSettings.fftWindow + 1)/ippgSettings.fftShiftSize);
  
  iPPG = zeros(ippgSettings.nMethod, nFrames);  
  hrEstimatedFFT = zeros(ippgSettings.nMethod, nAverageHRValues);
  spectrum = cell(1, ippgSettings.nMethod);  
  hrEstimatedWavelet = zeros(ippgSettings.nMethod, nAverageHRValues);

  settings = ippgSettings;
  for i = 1:ippgSettings.nMethod
    % compute iPPG from color signals
    settings.extractionMethod = ippgSettings.extractionMethod(i); 
    settings.processing = ippgSettings.processing(i);
    iPPG(i, :) = compute_ippg(channelsIntensity, settings);
    
    % estimate pulse rate for iPPG computed by all methods
    [hrEstimatedFFT(i, :), spectrum{i}] = DFT_pulse_rate_estimate(iPPG(i, :), settings);

    % we compute wavelet pulse rate estimates but do not evaluate them
    hrEstimatedWavelet(i, :) = wavelet_pulse_rate_estimate(iPPG(i, :), settings);
  end
  hrEstimated{iFile} = hrEstimatedFFT;
  
  % averaged reference pulse rates
  nWin = 1 + floor((nFrames - ippgSettings.fftWindow)/ippgSettings.fftShiftSize);
  hrTrue{iFile} = zeros(1, nWin);
  for iWin = 1:nWin       
    pos1 = 2*(iWin - 1) + 1;
    pos2 = 2*(iWin + 1);
    hrTrue{iFile}(iWin) = 1./mean(1./REFERENCE_PULSE_RATE{iFile}(pos1:pos2));  
  end

  if (RGB_DATA_FLAG(iFile) == 1) %we compare overall results only for RGB recordings
    subjectHRtrue{SUBJECT_INDEX(iFile)} = [subjectHRtrue{SUBJECT_INDEX(iFile)}, hrTrue{iFile}]; 
    subjectHRestimate{SUBJECT_INDEX(iFile)} = [subjectHRestimate{SUBJECT_INDEX(iFile)}, hrEstimated{iFile}(mainMethodIndex, :)];
    
    currentSectionIndex = SECTION_INDEX_PER_SUBJECT(iFile)*ones(1, length(hrTrue{iFile}));
    subjectSessionIndex{SUBJECT_INDEX(iFile)} = [subjectSessionIndex{SUBJECT_INDEX(iFile)}, currentSectionIndex];
  end
  
  % compute quality metrics for pulse rate estimation
  [shareErrorBelow3p5BPM{iFile}, shareErrorBelow7BPM{iFile}, ...
   corrCoef{iFile}, meanError{iFile}, ...
   rmse{iFile}, stdError{iFile}, corrPvalue{iFile}] = assess_estimation_performance(hrTrue{iFile}, hrEstimated{iFile});
    
  %compute SNR values and save only those for main method without wavelets
  dFreq = ippgSettings.samplingRate/ippgSettings.fftWindow;
  freqBPM = 60*(ippgSettings.minFreq:dFreq:ippgSettings.maxFreq);
  snrValues = compute_SNR(spectrum, freqBPM, hrTrue{iFile}, nBins{iFile});
  snr{iFile} = snrValues(:, noWaveletIndex);  % SNR is computed for non-wavelet filtered ppg
  
  xt{iFile} = ((ippgSettings.fftWindow/2)/ippgSettings.samplingRate)*(1:nWin);
  if (RGB_DATA_FLAG(iFile) == 1)
    fileID = fopen(['dataset' filesep filename{iFile} filesep filename{iFile} '_interframeDif.txt'],'r');
    motionData{iFile} = fscanf(fileID, '%f\n');
    fclose(fileID);

    errorForMotion{iFile} = -abs(hrEstimated{iFile}(mainMethodIndex, :) - hrTrue{iFile});
    startPosForMotion{iFile} = 1:ippgSettings.fftShiftSize:(length(motionData{iFile})-ippgSettings.fftWindow+1);
    endPosForMotion{iFile} = startPosForMotion{iFile}+ippgSettings.fftWindow;
    dFreqMotion{iFile} = dFreq;
  end 
end

%% compute percentage of epochs within 5% or 10% of subject mean pulse rate 
maxNSectionPerSubject = max(SECTION_INDEX_PER_SUBJECT);
shareErrorBelow5percent = zeros(1, nFile);
shareErrorBelow10percent = zeros(1, nFile);
nBelow5Percent = 0;
nBelow10Percent = 0;
for iSubject = 1:nSubject
  diffHR = abs(subjectHRestimate{iSubject} - subjectHRtrue{iSubject});
  nBelow5Percent = nBelow5Percent + nnz(diffHR < 0.05*subjectHRtrue{iSubject});
  nBelow10Percent = nBelow10Percent + nnz(diffHR < 0.1*subjectHRtrue{iSubject});
  for iSession = 1:maxNSectionPerSubject
    iFile = find((SUBJECT_INDEX == iSubject) & ...
                 (SECTION_INDEX_PER_SUBJECT == iSession));
    trialIndex = (subjectSessionIndex{iSubject} == iSession);
    if (nnz(trialIndex) > 0)  
      shareErrorBelow5percent(iFile) = nnz(diffHR(trialIndex) < 0.05*subjectHRtrue{iSubject}(trialIndex))/nnz(trialIndex);
      shareErrorBelow10percent(iFile) = nnz(diffHR(trialIndex) < 0.10*subjectHRtrue{iSubject}(trialIndex))/nnz(trialIndex);
    end  
  end
end
nTrial = sum(cellfun(@length, subjectHRtrue));
overallShareErrorBelow5percent = nBelow5Percent/nTrial;
overallShareErrorBelow10percent = nBelow10Percent/nTrial;
%disp(['Subject ' num2str(iSubject) ': ' num2str(length(subjectHRestimate{iSubject})) ' epochs, ' num2str(nCorrect) ' correct, ' num2str(nAccept - nCorrect) ' acceptable']);

%% report results
clc
disp(['Overall results: epochs with absolute error <5%: ' sprintf('%0.2f',overallShareErrorBelow5percent) ', <10%: ' sprintf('%0.2f',overallShareErrorBelow10percent) '.']);
   
for iFile = 1:nFile
  if (RGB_DATA_FLAG(iFile))  
    mainMethodIndex = 4;
  else
    mainMethodIndex = 1;
  end
  iSubject = SUBJECT_INDEX(iFile);
  iSessionPerSubject = SECTION_INDEX_PER_SUBJECT(iFile);
  nEpochs = length(hrEstimated{iFile});
  disp(['Session ' num2str(iFile) '(Subject ' num2str(iSubject) '): ' num2str(nEpochs) ' epochs, ']);
  disp(['  Mean absoulte error: ' sprintf('%0.2f',meanError{iFile}(mainMethodIndex)) ', Correlation: ' sprintf('%0.2f',corrCoef{iFile}(mainMethodIndex)) ';']);
  disp(['  Epochs with absolute error <3.5 BPM: ' sprintf('%0.2f',shareErrorBelow3p5BPM{iFile}(mainMethodIndex)) ', <7 BPM: ' sprintf('%0.2f',shareErrorBelow7BPM{iFile}(mainMethodIndex)) ';']);
  if (RGB_DATA_FLAG(iFile))  
    disp(['  Epochs with absolute error <5%: ' sprintf('%0.2f',shareErrorBelow5percent(iFile)) ', <10%: ' sprintf('%0.2f',shareErrorBelow10percent(iFile)) '.']);
  end  
end  

%% draw Bland Altman joint (Fig 5 from the manuscript)
figure
set( axes,'fontsize', FontSize, 'FontName', 'Times New Roman');

allHRestimate = [subjectHRestimate{1:nSubject}];
allHRtrue = [subjectHRtrue{1:nSubject}];
allSessionIndex = [];
addedSessionNumber = 0;
legendEntries = cell(1, 1);
for iSubject = 1:nSubject
  allSessionIndex = [allSessionIndex, subjectSessionIndex{iSubject} + addedSessionNumber];
  nSession = max(unique(subjectSessionIndex{iSubject}));
  for iSession = addedSessionNumber + 1:addedSessionNumber + nSession;
    legendEntries{iSession} = ['Session ' num2str(iSession)];
  end 
  addedSessionNumber = addedSessionNumber + nSession;
end  
disp('Bland-Altman plot: regression parameters.')
bland_altman_plot(allHRestimate, allHRtrue, allSessionIndex);
 
legend_handleMain = legend(legendEntries(:), 'location', 'NorthWest');
set(legend_handleMain, 'fontsize', FontSize, 'FontName', 'Times New Roman'); 
set( gca, 'fontsize', FontSize, 'FontName', 'Times New Roman');
xlabel( '(videoPR + refPR)/2 [BPM]', 'fontsize', FontSize, 'FontName', 'Times New Roman');
ylabel( 'videoPR - refPR [BPM]', 'fontsize', FontSize, 'FontName', 'Times New Roman');


%% motion plot (Fig 8 from the manuscript)
markerSize = 3;
motionFile = [3, 5];

imageLabel = {'A','B'};
yTickPlotwise = {[0.0020, 0.0021], 0.0025:0.0001:0.0030};
figure
set( axes,'fontsize', FontSize, 'FontName', 'Times New Roman');
for iPlot = 1:2
  iFile = motionFile(iPlot);
  
  nAverageMotionValues = min(length(xt{iFile}), length(startPosForMotion{iFile}));
  averageMotion = zeros(1, nAverageMotionValues);
  for i = 1:nAverageMotionValues
    averageMotion(i) = mean(motionData{iFile}(startPosForMotion{iFile}(i):endPosForMotion{iFile}(i)));
  end
  
  subplot(1,2,iPlot)
  yyaxis left
  plot(xt{iFile}, snr{iFile,1}, 'linewidth', LineWidth, 'marker', '*', 'markersize', markerSize);
  yyaxis right
  plot(xt{iFile}, averageMotion, 'linewidth', LineWidth, 'LineStyle', '--', 'marker', 'o', 'markersize', markerSize);
  axis tight;
  set( gca, 'fontsize', FontSize, 'FontName', 'Times New Roman');
  title(imageLabel{iPlot}, 'fontsize', FontSize, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
  xlabel( 'Time [s]', 'fontsize', FontSize, 'FontName', 'Times New Roman');
  yyaxis left
  ylabel('iPPG signal-to-noise ratio, [dB]', 'fontsize', FontSize, 'FontName', 'Times New Roman');  
  yyaxis right
  ylabel('Averaged interframe difference', 'fontsize', FontSize, 'FontName', 'Times New Roman', 'Interpreter', 'Tex');      
  set(gca, 'ydir', 'reverse', 'YTick', yTickPlotwise{iPlot}, 'YTickLabel', cellstr(num2str(yTickPlotwise{iPlot}(:))), 'YTickLabelRotation', 90);   
end  



