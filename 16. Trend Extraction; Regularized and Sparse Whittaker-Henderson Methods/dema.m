% dema.m - steady-state double exponential moving average 
%
% Usage: [a,b,a1,a2,cinit] = dema(y,N,cinit);
%        [a,b,a1,a2,cinit] = dema(y,N);        cinit = [y(1); 0], default
%
% y = signal to be smoothed (length L)
% N = equivalent length, N = (1+la)/(1-la) ==> la = (N-1)/(N+1), must have N<L
% cinit = any 2x1 vector of initial values of [a; b], i.e., the values at n=-1
%       = [y(1); 0], default, matches SEMA and corresponds to [a1; a2] = [y(1); y(1)] at n=-1
%       = 'f', fits straight line in first N samples, Refs.[1,2], cinit = [ones(N,1),(1:N)']\y(1:N)
%       = 'c', uses the value of the first EMA at time n=N-1 to serve as the initial value
%              for the second EMA starting at n=N, Ref.[3]
%
% a  = local level,         a = 2*a1-a2
% b  = local slope,         b = al/la*(a1-a2)
% a1 = output of 1st EMA,  a1 = a - la/al * b
% a2 = output of 2nd EMA,  a2 = a - 2*la/al * b
% cinit = numerical vector of initial conditions for a,b at n=-1
%
% notes: special case of MEMA with d=1,
%        implements explicitly the d=1 state-space STEMA version,
%        transient initialization effects disappear after about 6*N samples,
%        tau-step ahead prediction = a + b*tau,
%        for standard convolutional output, choose cinit = [0;0];
%        N is the equivalent SMA length, i.e., EMA and SMA have the same lag and NRR
%
% References: [1] D. C. Montgomery & L. A. Johnson, "Forecasting and Time Series Analysis", McGraw-Hill, NY, 1976.
%             [2] P. G. Mulloy, "Smoothing Data with Less Lag", Tech. Anal. Stocks & Commod., vol.12, p.11, Jan.1994. 
%             [3] S. B. Achelis, "Technical Analysis from A to Z", 2/e, McGraw-Hill, NY, 2001. 
%
% Example from [3], p.121, generated with 
%          y = [122.9060  126.5000  140.4060  174.0000  159.8120  170.0000  176.7500 ...
%               175.5310  166.5620  163.7500  170.5000  175.0000  184.7500  202.7810]';
%          N=5; [a,b,a1,a2,cinit] = dema(y,N,'c');  % emat(1,la)*cinit --> [a1in,a2in] at n=-1
%
%     n     y        a1        a2         a         b
%   -----------------------------------------------------
%    -1            122.9060  213.0267   32.7853  -45.0604   initial values
%     0  122.9060  122.9060  182.9865   62.8255  -30.0402
%     1  126.5000  124.1040  163.3590   84.8490  -19.6275
%     2  140.4060  129.5380  152.0853  106.9907  -11.2737
%     3  174.0000  144.3587  149.5098  139.2076   -2.5756 
%    *4  159.8120 *149.5098 *149.5098  149.5098    0.0000   * initialization of a2 at n=N-1
%     5  170.0000  156.3399  151.7865  160.8932    2.2767
%     6  176.7500  163.1432  155.5721  170.7144    3.7856
%     7  175.5310  167.2725  159.4722  175.0728    3.9001
%     8  166.5620  167.0357  161.9934  172.0780    2.5212
%     9  163.7500  165.9404  163.3090  168.5718    1.3157
%    10  170.5000  167.4603  164.6928  170.2278    1.3837
%    11  175.0000  169.9735  166.4530  173.4940    1.7602
%    12  184.7500  174.8990  169.2684  180.5297    2.8153
%    13  202.7810  184.1930  174.2432  194.1428    4.9749
% 
% see also mema, sema, tema, stema, wema, emat, lpbasis

% S. J. Orfanidis - 2009-2018

function [a,b,a1,a2,cinit] = dema(y,N,cinit)

if nargin==0, help dema; return; end

S = size(y); y = y(:);

if nargin==2, cinit = [y(1); 0]; end 

la = (N-1)/(N+1); al = 1-la;   % effective lambda
A = [1 1; 0 1];                % state matrix       
k = [1-la^2; (1-la)^2];        % Kalman gain vector
u1 = [1; 1];

if cinit == 'f'                        % fit line to first N data
   n = (1:N)';
   cinit = [n.^0, n] \ y(n);    
elseif cinit == 'c'                    % cascaded initialization
   a1in = y(1);                        % initialize first EMA 
   a1 = sema(y(1:N),N);                % calculate a1(N-1)
   a2in = a1(N);                       % initialize 2nd EMA at n=N-1
   for n = N:-1:1                      % iterate backwards to n=-1
      a2in = (a2in - al*a1(n))/la;     % determine a2in at n=-1
   end
%  alternative: a2=filter(-al/la,[1,-1/la],flip(a1),a1(N)/la); a2in=a2(end);
   cinit = emat(1,la) \ [a1in; a2in];  % transform to a,b basis
end                                    % note, emat(1,la) = [1,-la/al; 1,-2*la/al]
                                       
c1 = cinit;                    % initialize at n=-1
for n=1:length(y)              % run steady-state dema 
    e = y(n) - u1'*c1;         % prediction error
    c = A*c1 + k*e;            % steady-state dema
    c1 = c;                    % update state
    a(n) = c(1); 
    b(n) = c(2);
end

a1 = a - la/al * b;
a2 = a - 2*la/al * b;

% a1 = filter(al,[1,-la],y,la*a1in);     % alternative calculation, used in MEMA
% a2 = filter(al,[1,-la],a1,la*a2in);
% a = 2*a1-a2;
% b = al/la*(a1-a2);

a = reshape(a,S);
b = reshape(b,S);
a1 = reshape(a1,S);
a2 = reshape(a2,S);







