% acf.m - sample auto-correlation function
%
% Usage: Rxy = acf(x,x,M) = autocorrelation of length-N vector x
%        Rxy = acf(x,y,M) = cross-correlation between two length-N vectors x,y
%        Rxy = acf(x,y)   (equivalent to M=N-1)
%
% x = length-N vector 
% y = length-N vector
% M = maximum lag (typically, M <= N-1, pads zeros if M > N-1)
%
% Rxy = row vector of positive lags [Rxy(0), Rxy(1), ..., Rxy(M)]
%
% notes: calculates Rxy(k) = (1/N) \sum_{n=0}^{N-1-k} x(n+k) * conj(y(n)), k=0:M
%
%        x,y must have the same length N
%        x,y can be entered as row or column vectors
%        
%        use y=x for the "biased" sample autocorrelation of y
%
%        the double-sided cross-correlation Rxy(k), k=-M:M, can be calculated by
%        Rxy = acf(x,y,M); Ryx=acf(y,x,M); R = [conj(flipv(Ryx(2:end))), Rxy]
%        because Rxy(-k) = conj(Ryx(k))
%        
%        MATLAB's XCORR in R11.1 produces the complex-conjugate of Rxy(k)
%        but in version R13, it produces the same Rxy(k)
%
%        partial ACF, PARCORR, or reflection coefficients can be calculated by:
%
%             R = acf(y,y,M); gamma = lpg(lev(R));

% S. J. Orfanidis - 1999

function Rxy = acf(x,y,M)

if nargin==0, help acf; return; end

N = length(x);
x = x(:); 
y = y(:);

if nargin==2, M=N-1; end

Rxy = zeros(1,M+1);

if N==1, 
    Rxy(1) = x(1) * conj(y(1));
else
   for k=0:M,
       n = 0:N-1-k;
       Rxy(k+1) = x(n+k+1).' * conj(y(n+1)) / N;      
   end
end

