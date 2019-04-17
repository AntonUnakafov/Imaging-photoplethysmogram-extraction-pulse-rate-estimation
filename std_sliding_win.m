% @brief std_sliding_win computes running standard deviation of the signal 
% in windows of given length
% INPUT
%   - x - the signal, a row vector 1xN;
%   - w - the length of the sliding window. 
% OUTPUT:
%   - stdX - running standard deviation, a row vector 1xN
%     stdX(1)     is the std of x(1:w/2)...
%     stdX(w/2)   is the std of x(1:w)...
%     stdX(N-w/2) is the std of x(N-w+1:N)...
%     stdX(N)     is the std of x(N-w/2+1:N)...
%
% EXAMPLE OF USE:
%{
  %generate a composite signal of 200 points, first 100 having std 1, last 100 having std 2 
  signal = [randn(1, 100), 2*randn(1, 100)];
  w = 50;
  runningSTD = std_sliding_win(signal, w);
  figure
  plot(runningSTD)
  % the plot should have two plateus: y=1 for x = 1...75 and y=2 for x =
  125...200, with a transition in-between
%}

function stdX = std_sliding_win(x, w)
    % element count in each window
    n = conv(ones(1, length(x)), ones(1, w), 'same');
    
    s = conv(x, ones(1, w), 'same');
    q = x.^ 2;    
    q = conv(q, ones(1, w), 'same');
    stdX = sqrt(abs((q - s.^2 ./ n) ./ (n - 1))); %abs is used just in case
end