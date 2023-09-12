function [ysurr] = multivariate_surrogates_fix(y,n_surr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   input argument y where y(:,j) is the jth variable
%   input argument n_surr is the number of surrogates
%   output argument ysurr(:,j,k) is kth surrogates of jth variable
[L nvar] = size(y);
Y = fftshift(fft(y));

% columns now contain fft of each variable
% check for even or odd data
if(rem(L,2)==0)
    % even
    phase = 2*pi*rand(1,n_surr*(L/2-1)+1000);
    phase = phase(phase~=0);  %% don't double count 0 and 2*pi
    phase = phase(1:n_surr*(L/2-1));  %% reduce the amount of data
    phase = reshape(phase,L/2-1,n_surr);
    phase = [zeros(1,n_surr); phase; zeros(1,n_surr); -flipud(phase)];
    phase = exp(i*phase);
else 
    %odd
    phase = 2*pi*rand(1,n_surr*(L-1)/2+1000);
    phase = phase(phase~=0);  %% don't double count 0 and 2*pi
    phase = phase(1:n_surr*(L-1)/2);  %% reduce the amount of data
    phase = reshape(phase,(L-1)/2,n_surr);
    phase = [phase; zeros(1,n_surr); -flipud(phase)];
    phase = exp(i*phase);
end


for j=1:n_surr
    Ynew = Y.*phase(:,j);
    ysurr(:,:,j) = ifft(ifftshift(Ynew));
end   
