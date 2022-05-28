% sema.m - single exponential moving average 
%
% Usage: a = sema(y,N,ainit);
%        a = sema(y,N);          uses ainit = y(1)
%
% y = signal to be smoothed (length L)
% N = equivalent length, N = (1+la)/(1-la) ==> la = (N-1)/(N+1) 
% ainit = any initial value, defaults to ainit = y(1)
%
% a = smoothed version of y, represents local-level
%
% notes: single steady-state EMA with parameter la = (N-1)/(N+1)
%        N is the equivalent SMA length, 
%        i.e., EMA and SMA have the same lag and NRR
%
%        default initilization, ainit = y(1)
%        another popular choice is, ainit = mean(y(1:N))
%        for standard convolutional output, ainit = 0
%
% see also tech, sma, dema, tema, mema, stema, wema, gdema, t3, hma, ehma

% S. J. Orfanidis - 2009-2018

function a = sema(y,N,ainit)

if nargin==0, help sema; return; end
if nargin==2, ainit = y(1); end

la = (N-1)/(N+1);
 
aold = ainit;
for n=1:length(y)
    a(n) = la*aold + (1-la)*y(n);
    aold = a(n);
end

a = reshape(a,size(y));

