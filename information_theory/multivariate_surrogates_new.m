function [ysurr] = multivariate_surrogates(y,n_surr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   input argument y where y(:,j) is the jth variable
%   input argument n_surr is the number of surrogates
%   output argument ysurr(:,j,k) is kth surrogates of jth variable
[L nvar] = size(y);
Y = fftshift(fft(y));
% columns now contain fft of each variable
phase = 2*pi*rand(1,n_surr*L+1000);
phase = phase(phase~=0);
phase = phase(1:n_surr*L);
phase = reshape(phase,L,n_surr);
phase = exp(i*phase);
%phasearr = repmat(phase,1,nvar,1)
for j=1:n_surr
   Ynew = Y.*phase(:,j);
   Ynew = (Ynew + flipud(conj(Ynew)))/2;
   ysurr(:,:,j) = ifft(ifftshift(Ynew));
end   

