# Imaging-photoplethysmogram-extraction-pulse-rate-estimation
Tools for imaging photoplythysmogram extraction and processing

## Intro 
Imaging photoplethysmogram (iPPG) is a technique for remote non-contact pulse rate measurement. iPPG is usually acquired from facial or palm video.
This package provides tools for iPPG signal extraction and processing. The rhesus monkeys iPPG data from [1] are used as a test dataset.

**INPUT:** video file.

**OUTPUT:** iPPG signal; estimated pulse rate.

## Content
1. extract_color_channels_from_video extracts colour signal from the video. Colour signals are computed as values of red, green and blue colour component averaged over Region of Interest (ROI) for each video frame. ROI can be selected either manually for the first frame (if only limited amount of motion is expected) or set automatically using the Viola-Jones algorithm (for iPPG extraction from human face only!). This function optionally excludes from ROI non-skin and corrupted pixels.

2. compute_ippg implements iPPG extraction methods considered in [2] (including the recently introduced CHROM and POS methods) as well as some iPPG pre- and post-processing techniques.

3. ippg_extraction_example- basic (minimal) example of using the package for estimation of pulse rate from the iPPG extracted from a video.

4. dataset_analysis - extended example of using the package for the data from [1].

5. Signal processing techniques implemented as separate m-files: wavelet_filter, wavelet_init_scales, smoothness_priors_detrending, std_sliding_win.

6. Functions for pulse rate estimation from iPPG signal:

6.1. DFT_pulse_rate_estimate uses Discrete Fourier Transform to compute average pulse rate.

6.2. wavelet_pulse_rate_estimate uses Continuous Wavelet Transform to estimate pulse rate.

7. Useful functions for comparing iPPG-based pulse rate with the ground truth:

7.1. bland_altman_plot - draws Bland Alman plot for the data.

7.2. compute_SNR - computes Signal-to-Noise Ratio (SNR) of the iPPG signal given the true pulse rate.

7.3. assess_estimation_performance - computes a number of estimation quality metrics, including root-mean-square error, mean absolute error, Pearson correlation, etc.

8. Dataset folder contains the dataset used for testing the package. The dataset was recorded from rhesus monkeys, therefore pulse rate is higher than for humans (100-250 BPM), please refer to [1] for details.

9. dataset_description.docx contains a brief description of the dataset.

Additional functionality will be added later.

## Acknowledgements
I would like to thank Gasper Slapničar for extensive testing and improving the code.

I would like to thank Dr. Cardoso for posting [the jadeR.m script](http://www.mikexcohen.com/lecturelets/eigen/jader.m), which I use to implement ICA-based iPPG extraction.

## References
[1] [Unakafov AM, Moeller S, Kagan I, Gail A, Treue S, Wolf F. Using imaging photoplethysmography for heart rate estimation in non-human primates. PLoS ONE 2018;13(8): e0202581.](https://doi.org/10.1371/journal.pone.0202581)

[2] [Unakafov AM. Pulse rate estimation using imaging photoplethysmography: generic framework and comparison of methods on a publicly available dataset. Biomedical Physics & Engineering Express. 2018;4(4):045001.](https://doi.org/10.1088%2F2057-1976%2Faabd09)

## Cite As

Unakafov AM. Imaging photoplethysmogram extraction&pulse rate estimation https://www.mathworks.com/matlabcentral/fileexchange/67527, MATLAB Central File Exchange (2018). Retrieved December 12, 2018.

Unakafov AM. Pulse rate estimation using imaging photoplethysmography: generic framework and comparison of methods on a publicly available dataset. Biomedical Physics & Engineering Express. 2018;4(4):045001. https://doi.org/10.1088%2F2057-1976%2Faabd09

Unakafov AM, Moeller S, Kagan I, Gail A, Treue S, Wolf F. Using imaging photoplethysmography for heart rate estimation in non-human primates. PLoS ONE 2018;13(8): e0202581. https://doi.org/10.1371/journal.pone.0202581
