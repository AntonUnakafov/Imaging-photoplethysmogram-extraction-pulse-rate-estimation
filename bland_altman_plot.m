% @brief bland_altman_plot plots a Bland-Altman Plot, characterizing similarity
% of two estimates.
%
% INPUT
%   - estimate1, estimate2 - two estimates of a certain quantity in several 
%     epochs, row or column vectors of equal size;
%
% OPTIONAL INPUT
%   - sessionIndex - array of session indices of the same size as estimates,
%     for each epoch specifies to which session it belongs. 
%     This allows to plot epochs from different sessions in different colors       
%
% OUTPUT:
%   - means - the means of the estimates
%   - diffs - the differences of the estimates
%   - meanDiff - the mean difference (bias) of the estimates
%   - confInterval - the limits of the 95% confidence interval
%   - linRegressCoef - the coefficients of the linear regression between
%     means and diffs                 
function [means, diffs, meanDiff, confInterval, linRegressCoef] = bland_altman_plot(estimate1, estimate2, sessionIndex)
  colorList = [0.0, 0.0, 0.6; ...
               0.6, 0.5, 1.0; ...
               
               0.7, 0.4, 0.2; ...
               1.0, 0.8, 0.5; ...
               
               0.9, 0.0, 0.0; ...
               0.6, 0.1, 0.4; ...
               1.0, 0.6, 0.8; ...
               
               0.0, 0.9, 0.7; ...
               
               0.8, 1.0, 0.2; ...                   
               0.5, 1.0, 0.3; ...
               0.1, 1.0, 0.7; ...                 
               0.0, 1.0, 0.8; ...
               0.2, 0.6, 0.6; ...
               0.2, 1.0, 1.0; ...
               0.4, 1.0, 1.0; ...
               0.6, 0.9, 0.9; ...
               1.0, 0.2, 0.5; ...               
               0.7, 1.0, 1.0];
    
  if nargin==2
     sessionIndex = 1;
  end
    
  means = mean([estimate1;estimate2]);
  diffs = estimate1 - estimate2;
   
  % compute linear regression coefficients
  mdl = fitlm(means, diffs)
  linRegressCoef = polyfit(means, diffs, 1); 

  meanDiff = mean(diffs);
  sdDiff = std(diffs);
  confInterval = [meanDiff + 1.96*sdDiff, meanDiff - 1.96*sdDiff]; % 95% confidence interval
        
  minX = min(means) - 2;
  maxX = max(means) + 2;
  minY = min(min(diffs) - 2, confInterval(2) - 3);
  maxY = max(max(diffs) + 2, confInterval(2) + 3);
    
  % plot results
  hold on
  nColor = max(unique(sessionIndex));
  for iColor = 1:nColor
    index = (sessionIndex == iColor);
    plot(means(index), diffs(index), 'o', 'markersize', 4, 'Color', colorList(iColor, :));
  end
    
  % plot mean and confidence interval
  plot([minX maxX], [meanDiff meanDiff], 'k-');
  plot([minX maxX], [confInterval(1) confInterval(1)],'k--'); 
  plot([minX maxX], [confInterval(2) confInterval(2)],'k--'); 
 % text(minX + 2, confInterval(1) - 1, 'Mean + 1.96SD','fontsize', 9, 'FontName', 'Times');
 % text(minX + 2, confInterval(2) - 1, 'Mean - 1.96SD','fontsize', 9, 'FontName', 'Times');

 % plot the regression line
  plot([minX maxX], [minX maxX]*linRegressCoef(1)+linRegressCoef(2),'b:', 'linewidth', 1.5); 
  hold off
  axis ( [minX, maxX, minY - 1, maxY + 1 + 0.15*(maxY - minY)] );
 % axis ( [minX, maxX, -16.5, 28.5] );
 % axis ( [minX, maxX, -35.5, 34.5] );
end