% @brief extract_color_channels_from_video extracts the color channels
% from the color video. See [a] for details
%
% CITING THE CODE
% [a] Unakafov AM. Pulse rate estimation using imaging photoplethysmography: generic framework and comparison of methods on a publicly available dataset. Biomedical Physics & Engineering Express. 2018;4(4):045001.
%
% INPUT
%   - filename - name of video file
%   - videoSettings - structure with video processing settings. 
%     It must contain the following fields:
%       - cascadeROIselection - boolean, whether face detection by the 
%         Viola-Jones algorithm [1] is used. 
%         This requires Computer Vision toolbox!
%         Currently only as beta-version!!!
%       - manualROIselection - boolean, whether ROI is selected manually
%         for the first frame. Can be used in combination with previous
%         setting, then face is searched only in the specified ROI.
%         average pulse rate
%       - isSTDmaskingOn - boolean, specifies whether outliers elimination is applied
%         Outliers determined by comparison with std of a color channel
%         multimlied by a constant stdCoef, see [a] for details
%       - stdCoef - coefficient for outliers elimination, used when
%         isSTDmaskingOn == true
%       - isHSVmaskingOn - boolean, specifies whether skin pixels are selected
%         using thresholding of HSV values.
%       - hsvMin - vector of 3 double values in range 0-1, specifies minimal 
%         values of hue, saturation and value for non-skin pixels filtering
%       - hsvMax - vector of 3 double values in range 0-1, specifies maximal 
%         values of hue, saturation and value for non-skin pixels filtering
%
% OUTPUT:
%   - rawColorSignal - matrix 3xN containing three color channels (RGB) extracted from video
%
% REFERENCES: 
% [1] Viola P, Jones M (2001) 
% Proceedings of the 2001 IEEE Computer Society Conference on 
% Computer Vision and Pattern Recognition 1 I-511
%

function rawColorSignal = extract_color_channels_from_video(filename, videoSettings)
videoObject = VideoReader(filename);

% filling the missing fields of the settings stucture 
if (~isfield(videoSettings,'cascadeROIselection'))
    videoSettings.cascadeROIselection = false;
end
if (~isfield(videoSettings,'manualROIselection'))
    videoSettings.manualROIselection = false;
end

nFrame = ceil(videoObject.FrameRate*videoObject.Duration);

if (hasFrame(videoObject)) % set initial ROI
    currentFrame = readFrame(videoObject);
    if (videoSettings.manualROIselection)
        % manual selection of ROI rectangle in the first frame
        rectROI = loc_manual_ROI_selection(currentFrame);
    else
        % ROI rectangle covers the whole frame
        rectROI = [1,1,videoObject.Width, videoObject.Height];
    end 
    widthROI = rectROI(3);
    heightROI = rectROI(4);     
    rawColorSignal = zeros(3, nFrame); % preallocate array
else
    error('Video processing failed: no frames are available!');
end

if (videoSettings.cascadeROIselection) % initialize face detector
    faceDetector = vision.CascadeObjectDetector;
    faceExpectedRect = [1, 1, widthROI, heightROI];
end    

% process all individual frames
for iFrame = 1:nFrame
    % disp(strcat(['Processing frame:', ' ', num2str(iFrame)]))
    if (iFrame > 1) %the first frame is already loaded
        if (~hasFrame(videoObject))
            break
        end   
        currentFrame = readFrame(videoObject);
    end
    croppedImage = loc_crop_image(currentFrame, rectROI); % crop initial ROI rectangle

    if (videoSettings.cascadeROIselection) % detect face if necessary
        faceRect = loc_cascade_ROI_selection(croppedImage, faceExpectedRect);     
        facialImage = loc_crop_image(croppedImage, faceRect);
    else
        facialImage = croppedImage;        
    end
    currentROI = loc_ROI_refinement(facialImage, videoSettings); % filter out bad pixels
    
    % check whether ROI due to face detection didn't become less than 4 pixels
    if (videoSettings.cascadeROIselection) 
        if (sum(sum(currentROI)) < 4) 
            facialImage = loc_crop_image(croppedImage, faceExpectedRect);
            currentROI = loc_ROI_refinement(facialImage, videoSettings);
            faceRect = faceExpectedRect;
        else
            faceExpectedRect = faceRect;
        end    
    end    
    
    % show the video frames
    maskedImage = facialImage;
    for iChannel = 1:3
        colorChannel = maskedImage(:,:,iChannel);
        colorChannel(~currentROI) = 1;
        maskedImage(:,:,iChannel) = colorChannel;
    end     
    imshow(maskedImage);
    title(['Time: ' num2str(iFrame/videoObject.FrameRate) ' / ' num2str(videoObject.Duration)]);
    drawnow update
    
    % compute the average color values over the ROI 
    for iChannel = 1:3
        colorChannel = facialImage(:,:,iChannel);
        rawColorSignal(iChannel,iFrame) =  mean(mean(colorChannel(currentROI)));
    end
end
end

% manual selection of ROI rectangle in the first frame
function rectROI = loc_manual_ROI_selection(currentFrame)
figureHandle = figure('Name', 'Selecting the region of interest');
imshow(currentFrame);
title('Please select the region of interest with a mouse:');
h = imrect;
if isempty(h)
    warning('Video processing interrupted by user');
    return;
end
h.Deletable = false;
addNewPositionCallback(h,@(p) title([mat2str(round(p),3), ', double click on ROI to continue...']));
fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
setPositionConstraintFcn(h,fcn);
title([mat2str(round(getPosition(h)),3), ', double click on ROI to continue...']);
rectROI = round(wait(h));
if isempty(rectROI)
    warning('Video processing interrupted by user');
    return;
end
close(figureHandle);
end

% localising the face using Viola-Jones algorithm
function faceRect = loc_cascade_ROI_selection(currentFrame, faceExpectedRect)

faceDetector = vision.CascadeObjectDetector;
[~, frameWidth, ~] = size(currentFrame);

faceRect = faceExpectedRect; % if we do not find a new face, simply use the last location

faces = step(faceDetector, currentFrame);
[nFace, ~] = size(faces);

if (nFace > 0) % faces found
    % center of faceRect cannot move between successive frames on more than a half of the faceRect diagonal
    maxDist = sqrt(faceExpectedRect(3)^2 + faceExpectedRect(4)^2)/2;
    % area cannot decrease more than twice
    minArea = faceExpectedRect(3)*faceExpectedRect(4)/2;
    if (faceExpectedRect(3) >= 0.95*frameWidth) % if expected face is too large
        minArea = 9;  % consider it incorrect and set facial area threshold to a small const
    end
    expectedCenter = [faceExpectedRect(1) + 0.5*faceExpectedRect(3), ...
        faceExpectedRect(2) + 0.5*faceExpectedRect(4)];
    
    % check all detected faces to find the one 
    for iFace = 1:nFace
        faceCenter = [faces(iFace, 1) + 0.5*faces(iFace, 3), ...
            faces(iFace, 2) + 0.5*faces(iFace, 4)];
        currDist = norm(faceCenter - expectedCenter);
        if ((currDist < maxDist) && (faces(iFace, 3)*faces(iFace, 4) > minArea))
            maxDist = currDist;
            faceRect = faces(iFace, :);
            % reduce the width of the face as it is done by Picard group
            faceRect(1) = faceRect(1) + 0.1*faceRect(3);
            faceRect(3) = 0.8*faceRect(3);
        end
    end
end


end      


% function for filtering bad pixels (non-skin or corrupted with artefacts)
function currentROI = loc_ROI_refinement(croppedImage, videoSettings)
    [heightROI,widthROI,~] = size(croppedImage);
    currentROI = ones(heightROI,widthROI);
    % HSV filtering: pixels with hue, saturation or value 
    % outside of the specified range are discarded
    if (videoSettings.isHSVmaskingOn && (length(videoSettings.hsvMin) == 3) && (length(videoSettings.hsvMax) == 3))
        hsvFrame = rgb2hsv(croppedImage);        
        for iChannel = 1:3
            hsvMask = roicolor(hsvFrame(:,:,iChannel), videoSettings.hsvMin(iChannel), videoSettings.hsvMax(iChannel));    
            currentROI = currentROI & hsvMask;
        end            
    end    
    
    % std thresholding: pixels with color channels too far from the mean
    % are discarded
    if (videoSettings.isSTDmaskingOn)
        stdMask = currentROI;
        for iChannel = 1:3
            colorChannel = croppedImage(:,:,iChannel);
            channelMean = mean(mean(colorChannel(currentROI)));
            channelStd  = std2(colorChannel(currentROI));
            minChannelValue = channelMean - videoSettings.stdCoef*channelStd; 
            maxChannelValue = channelMean + videoSettings.stdCoef*channelStd;
            % we do not modify here currentROI since we need it to compute std corretly
            stdMask = stdMask & roicolor(colorChannel, minChannelValue, maxChannelValue);                
        end 
        currentROI = currentROI & stdMask;
    end
end

% wrapper for cropping the image
function croppedImage = loc_crop_image(image, rect)
    croppedImage = imcrop(image, [rect(1:2), rect(3:4)-1]);        
end        