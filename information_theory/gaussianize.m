function [yg ysorted] = gaussianize(y)
%UNTITLED2 Summary of this function goes here
%   y is a multivariate time series with each
%   variable's time series in each column
%   yg returns a Gaussianized time series
%   b contains the sorted variables
[L nvar] = size(y);
[ysorted ind] = sort(y);
[~, ind2] = sort(ind);
g = randn(1,L);
[bg gind] = sort(g);
yg = bg(ind2);
end

