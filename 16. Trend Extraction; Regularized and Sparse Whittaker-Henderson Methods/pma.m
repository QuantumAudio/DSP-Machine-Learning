% pma.m - predictive moving average, linear fit
%
% Usage: y = pma(x,N,tau,yin)
%
% x = signal to be filtered
% N = filter length, filter order N-1
% tau = prediction distance, need not be an integer, (tau>0 or tau<0 corresponds to prediction or delay) 
% yin = 'f', progressive filtering for first N-1 outputs, default
%     = 'n', first N-1 outputs are NaNs, simplifies plotting
%     = 'c', standard convolutional output transients
%
% y  = predicted output, i.e., y(n) is the prediction of x(n+p) based on x(n),x(n-1),...,x(n-N+1)
%
% notes: predictive causal FIR filters h(n) are the minimum-norm solution of the constraint equations:
%           for d=1, \sum_{n=0}^{N-1} [1, n] * h(n) = [1, -tau] 
%
%        filter h has prescribed filter delay = nbar = \sum_{n=0}^{N-1} n*h(n) = -tau
%
%        some special cases equivalent to some Technical Analysis tools are, with d=1:
%           tau = -(N-1)/2,   ==> SMA, filter delay = (N-1)/2
%           tau = -(N-1)/3,   ==> WMA, filter delay = (N-1)/3
%           tau = 0,          ==> "Linear Regression Indicator", filter delay = 0
%           tau = 1,          ==> "Time Series Forecast", filter delay =-1, also used in the "Forecast Oscillator"
%           (tau=1)-(tau=0),  ==> "Linear Regression Slope"
%           Tech. Anal. Reference: S. B. Achelis, "Technical Analysis from A to Z", 2/e, McGraw-Hill, NY, 2001. 
%
% examples:      x = [10    12    11    13    16        17        19        17]
%       pma(x,5,1) = [NaN   NaN   NaN   NaN   16.3000   18.3000   21.2000   19.7000]
%       pma(x,5,0) = [NaN   NaN   NaN   NaN   15.0000   16.8000   19.2000   18.6000]
%    pma(x,5,-4/2) = [NaN   NaN   NaN   NaN   12.4000   13.8000   15.2000   16.4000] = sma(x,5)
%    pma(x,5,-4/3) = [NaN   NaN   NaN   NaN   13.2667   14.8000   16.5333   17.1333] = wma(x,5)

% S. J. Orfanidis - 2009-2018

function y = pma(x,N,tau,yin)

if nargin==0, help pma; return; end
if nargin<=3, yin = 'f'; end             % default 

S = size(x); x = x(:);

h = pmaimp(N,tau);           % PMA impulse response 

y = filter(h,1,x);        

if yin=='f'                  % filters of increasing order
    y(1) = x(1);             % no filtering if N=1
    for K=2:N-1
        hk = pmaimp(K,tau);
        yk = filter(hk,1,x(1:K));
        y(K) = yk(end);
    end
end

if yin=='n'                  % N-1 outputs are NaNs
    y(1:N-1) = NaN(N-1,1);
end

% else, y(1:N-1) = standard convolutional transients
             
y = reshape(y,S);

% --------------------------

function h = pmaimp(N,tau)

n=0:N-1;

h = 2*(2*N-1-3*n)/N/(N+1) + tau * 6*(N-1-2*n)/N/(N^2-1);

if N==1, h=1; end           % avoids NaN if N=1






