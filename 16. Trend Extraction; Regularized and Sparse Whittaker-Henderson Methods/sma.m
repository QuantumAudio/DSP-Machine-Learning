% sma.m - simple moving average
%
% Usage: y = sma(x,N,yin)
%
% x = signal to be smoothed
% N = filter length, filter order M=N-1
%
% yin = 'f', progressive filtering for first N-1 outputs, default
%     = 'n', first N-1 outputs are NaNs, simplifies plotting
%     = 'c', standard convolutional output transients
%
% y  = smoothed output 
%
% notes: plain averaging filter h(n) = 1/N, n=0,1,...,N-1    

% S. J. Orfanidis - 2009-2018

function y = sma(x,N,yin)

if nargin==0, help sma; return; end
if nargin==2, yin = 'f'; end             % default 

S = size(x); x = x(:);

h = ones(N,1)/N;

y = filter(h,1,x);           % standard convolutional output 

if yin=='f'                  % filters of increasing order
    y(1) = x(1);             % no filtering if N=1
    for K=2:N-1
        hk = ones(K,1)/K;
        yk = filter(hk,1,x(1:K));
        y(K) = yk(end);
    end
end

if yin=='n'                  % N-1 outputs are NaNs
    y(1:N-1) = NaN(N-1,1);
end

% else, y(1:N-1) = standard convolutional transients
           
y = reshape(y,S);


