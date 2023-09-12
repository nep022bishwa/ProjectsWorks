function [wM,bV] = caaft(xV, M, tau_max, p, K)
% [wM,bV] = caaft(xV, M, tau_max, p, K)
% CAAFT (Corrected Amplitude Adjusted Fourier Transform )
% generates surrogate data for a given time series 'xV' using the CAAFT  
% algorithm (statically transformed realisations of a Gaussian
% (autoregressive) process, which have the same amplitude distribution 
% and autocorrelation as the given time series 'xV').
% The generated surrogate data are like AAFT, but corrected to match 
% the autocorrelation.
% For the correction, a linear interpolation for the graph of 
% the relation between the Gaussian and the transformed
% autocorrelation is found, for lags up to a given 'tau_max'. Using this 
% interpolation function, for the autocorrelation of the given 'xV',
% the autocorrelation of the respective Gaussian time series is 
% estimated. Based on this autocorrelation, the coefficients of the 
% corresponding AR model of order 'p' are estimated and an AR-time 
% series is generated, and transformed to match the amplitude 
% distribution of 'xV'. 
% This process is done 'K' times to obtain a statistic of
% candidate AR models. Then the most proper is selected, in the sense 
% that the autocorrelation of the generated surrogate matches best the 
% autocorrelation of 'xV'. Based on this AR model, 'M' realizations 
% are generated and transformed to match the amplitude distribution of 
% 'xV'.
% 
% INPUT:
% - xV      : The original time series (column vector)
% - M       : The number of surrogate data to be generated 
% - tau_max : The maximum lag for autocorrelation
% - p       : The order of the AR model
% - K       : The number of candidate AR models to be estimated
% OUTPUT:
% - wM      : The CAAFT surrogate data (matrix of 'M' columns)
% - bV      : The parameters of the core AR model (column vector) 
% 
% NOTE: 
% If the fifth argument is not specified, then 'K' = 'M'.
% If the fifth and fourth argument are not specified, then 
% 'K' = 'M' and 'p = tau_max'.
% The functions 'levinson' and 'polystab' of the Signal Processing
% Toolbox are called.

if nargin == 4
  K = M;
elseif nargin == 3
  K = M;
  p = tau_max;
end
if tau_max<p
  tau_max = p;
end
ntrans = 100;   % transients for the generation of data from AR model
ntries = 10;    % Number of tries to be made seeking for stable AR model 
n = length(xV);

% The following gives the rank order, ixV, and autocorrelation, autxV, of xV  
[oxV,T1] = sort(xV);
[T,ixV]=sort(T1);
autxV = xcorr(xV-mean(xV),tau_max,'coeff');
autxV(1:tau_max) = [];

% In the following loop, the algorithm is repeated 'K' times. 
% The AR model coefficients and the autocorrelation of the surrogates
% are stored in the respective matrices (bM and autwM).  
bM = zeros(p+1,K);
autwM = zeros(tau_max+1,K);
for sample=1:K
  disp([' Sample ',int2str(sample),' for CAAFT, p=',int2str(p),' tau_max=',int2str(tau_max)])

  successflag = 0; % being 1 if a stabilised AR is found
  tries = 1;  
  while tries<=ntries & ~successflag
    % ===== AAFT algorithm and autocorrelations
    rV = randn(n,1); % Rank ordering white noise with respect to xV 
    [orV,T]=sort(rV);
    yV = orV(ixV);
    % >>>>> Phase randomisation (Fourier Transform): yV -> yftV 
    if rem(n,2) == 0
      n2 = n/2;
    else
      n2 = (n-1)/2;
    end
    tmpV = fft(yV,2*n2);
    magnV = abs(tmpV);
    fiV = angle(tmpV);
    rfiV = rand(n2-1,1) * 2 * pi;
    nfiV = [0; rfiV; fiV(n2+1); -flipud(rfiV)];
    % New Fourier transformed data with only the phase changed
    tmpV = [magnV(1:n2+1)' flipud(magnV(2:n2))']';
    tmpV = tmpV .* exp(nfiV .* i); 
    % Transform back to time domain
    yftV=real(ifft(tmpV,n));
    % <<<<<
    [T,T2] = sort(yftV); % Rank ordering xV with respect to yftV 
    [T,iyftV] = sort(T2);
    zV = oxV(iyftV);  % zV is the AAFT surrogate of xV
    % Autocorrelation of the phase randomised rank ordered noise, yftV,
    % and the AAFT surrogate, zV.
    autyftV = xcorr(yftV-mean(yftV),tau_max,'coeff');
    autyftV(1:tau_max) = [];
    autzV = xcorr(zV-mean(zV),tau_max,'coeff');
    autzV(1:tau_max) = []; 

    % ===== Linear interpolation function for autocorrelations
    % Order autzV for interpolation (the x-input of function 'interp1'
    % must be at increasing order) and match the rank order of autyft
    % accordingly
    [oautzV,iautzV] = sort(autzV); 
    oautyftV = autyftV(iautzV); 
    autvV = interp1(oautzV,oautyftV,autxV,'linear');
    % For values of autxV < min(autzV), [note that min(autzV)=oautzV(1)]
    % the interpolation is not defined (matlab gives NaN). 
    % Extrapolate using the ratio autyft/autz, i.e. the simple formula: 
    %    autvV(i) = autxV(i) * oautyftV(1) / oautzV(1)
    extV = find(isnan(autvV)==1);
    for ii=1:length(extV)
      autvV(extV(ii)) = (oautyftV(1)/oautzV(1))*autxV(extV(ii));
    end

    % ===== Derivation of stabilised AR model from autocorrelation 
    % Use the Levinson recursive algorithm to find the AR parameters
    [bV,E]=levinson(autvV,p);
    % If needed, find new bV with all roots in the unit circle. 
    if any(abs(roots(bV)) >= 1)
      bV = polystab(bV);
      if all(abs(roots(bV)) < 1)
        successflag = 1; 
      end   
    else
      successflag = 1;
    end

    % ===== Generation of the transformed normal time series
    % If a stabilised AR model is found continue to generate AR data
    if successflag  
      bV=bV(:);
      uV = zeros(n+p+ntrans,1);
      uV(1:p) = randn(p,1);
      for ii=p+1:n+p+ntrans
        uV(ii) = - bV(2:p+1)' * uV(ii-1:-1:ii-p) + randn;
      end
      uV = uV(p+ntrans+1:n+p+ntrans);
      % Transform the normal time series to match the amplitude 
      % distribution of the original data
      [tmp1V,tmp2V] = sort(uV);
      [tmp1V,iuV] = sort(tmp2V);
      wV = oxV(iuV); 
    else
      tries = tries+1;
    end
  end
  if tries>ntries
    error(['Could not find stable AR(',int2str(p),') model.'])
  end

  bM(:,sample)=bV;
  autwV = xcorr(wV-mean(wV),tau_max,'coeff');
  autwV(1:tau_max) = [];
  autwM(:,sample) = autwV;
end

% ===== Selection of the "most appropriate" AR model
% The "most appropriate" in the sense that the autocorrelation of the
% respective trasnformed AR time series is "closest" to the original 
% of all other candidates. Here, "closest" is defined as follows:
% 1. For each lag tau=1...tau_max, assume the autocorrelation vectors of 
%    lag components (1,...,tau), one for each of the 'sample' candidates.
% 2. Find the vector (of length tau) closest to the original 
%    autocorrelation vector (for lags 1,...,tau) using the Euclidean norm.
% 3. Repeat step 2 for tau=1 to tau_max. Find the candidate, which succeeds 
%    minimum distance most of the 'tau_max' times.   
candV = zeros(tau_max,1);
for tau=1:tau_max
  tauV = [1:tau]'; % Note that lag=0 is not included
  dautwM = autxV(tauV+1) * ones(1,K) - autwM(tauV+1,:);
  disautV = zeros(K,1);
  for jj=1:K
    disautV(jj) = norm(dautwM(:,jj));
  end
  [dbest,candV(tau)]=min(disautV);     
end
ocandV = sort(candV);
indV = ocandV(1);   % Contains the indices for minima without repetition
countV = 1;         % Contains the number the respective index is repeated
inow = 1;
for ii=2:tau_max
  if ocandV(ii) == indV(inow)
    countV(inow) = countV(inow)+1;
  else
    inow = inow+1;
    indV = [indV; ocandV(ii)];
    countV = [countV; 1];
  end
end
[tmpV,ii] = max(countV);
jbest = indV(ii);

% ===== Generation of surrogates with the "best" AR-model:
bV = bM(:,jbest); 
wM = zeros(n,M);
for ii=1:M
  uV = zeros(n+p+ntrans,1);
  uV(1:p) = randn(p,1);
  for jj=p+1:n+p+ntrans
    uV(jj) = - bV(2:p+1)' * uV(jj-1:-1:jj-p) + randn;
  end
  uV = uV(p+ntrans+1:n+p+ntrans);
  [tmp1V,tmp2V] = sort(uV);
  [tmp1V,juV] = sort(tmp2V);
  wM(:,ii) = oxV(juV); 
end

