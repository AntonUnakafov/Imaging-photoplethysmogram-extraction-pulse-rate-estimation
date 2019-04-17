% @brief assess_estimation_performance reports how well provided estimates 
% fit to the reference (ground truth) data
%
% INPUT
%   - refValues  - vector 1xN of reference (ground truth) values
%   - estimatedValues - matrix MxN containing M different estimates for each value;
%
% OUTPUT:
%   - shareErrorBelow3p5 - vector 1xM, percentage of values estimated with error below 3.5
%     (for each of the M estimation methods)
%   - shareErrorBelow7 - vector 1xM, percentage of values estimated with error below 7
%   - corrCoef - vector 1xM, the Pearson correlation coefficient betweeen
%     estimated and reference values
%   - meanError - vector 1xM, mean absolute estimation error
%   - rmse - vector 1xM, root mean square estimation error
%   - stdError - standard deviation of the absolute estimation error
%   - corrPvalue - vector 1xM, P-values for the Pearson correlation coefficients
%
function [shareErrorBelow3p5, shareErrorBelow7, corrCoef, meanError, rmse, ...
          stdError, corrPvalue] = assess_estimation_performance (refValues, estimatedValues)
  [nMethods, N] = size(estimatedValues);    
  error = repmat(refValues, nMethods, 1); % nMethods copies of refValues  
  error = abs(error - estimatedValues);  % absolute error between two estimates
  shareErrorBelow3p5 = zeros(1, nMethods);
  shareErrorBelow7 = zeros(1, nMethods);
  corrCoef = zeros(1, nMethods);
  corrPvalue = zeros(1, nMethods);
  for i = 1:nMethods
    shareErrorBelow3p5(i) = length(find(error(i, :) < 3.5))/N;
    shareErrorBelow7(i) = length(find(error(i, :) < 7.0))/N;   
    [corrCoef(i), corrPvalue(i)] = corr(refValues', estimatedValues(i, :)');
  end
  meanError = mean(error, 2)';
  stdError = std(error, 0, 2)';
  rmse = sqrt(mean(error.*error, 2))'; 
end
