%% 1. Sparse Spike Deconvolution
load yspike.mat

%% 1.a.Input Signal, Noisy Observations, and Filter
N = length(y);
M = 53;
n0 = 25;
L  =  N - M;
n = (1:L)';
Ln = (0:L-1)';
Mn = (0:M)';
Yn = (1:N)';

h = cos(0.15*(Mn - n0)).*exp(-0.004*(Mn - n0).^2); % define our filter
Hn = (1:length(h))';

w = linspace(0,pi,1000);

HH = 20*log10(abs(freqz(h,1,w))); % calculate magnitude response

H = convmtx(h(:),L);    % construct convolution matrix

s = zeros(L,1); % construct sparse input for comparison to recovered signal

ni = [20, 40 , 60, 70 ,80, 100, 120, 140]';

ai = [10, 8, 4, -4, 5, 6, -2, 4]';

s(ni) = ai;

% Plots

figure(1)
subplot(2,2,1)
plot(n,s,'r'), xlabel('n'), title('exact input, s(n)'),
axis([0, 200, -5, 11])

subplot(2,2,2)
plot(Yn,y), xlabel('n'), title('noisy obserbations, y(n)'),
axis([0, 200, -5, 11]), text(140,9,'SNR = 38 dB')

subplot(2,2,3)
plot(Hn,h), xlabel('n'), title('impulse response h(n)'),
axis([0, 200, -0.4, 1.2]), grid on

subplot(2,2,4)
plot(w/pi,HH), xlabel('\omega/\pi'), ylabel('dB'),
title('magnitude response in dB, |H(\omega)|'),
axis([0, 1, -40, 40]), grid on

%% 1.b. Ordinary Least Squares and L2-Regularized Solutions, lambda = 0.01

% Ordinary Least Squares Solution
xord = (H'*H)\(H'*y);

PerrorXord = 100*norm(xord - s)/norm(s);

% L2-Regularized Solution
lambda = 0.01;
x1L2 = (lambda*speye(L) + H'*H)\H'*y;

PerrorXL2 = 100*norm(x1L2 - s)/norm(s);

% Plots

figure(2)
subplot(1,2,1)
plot(n,xord,'r'), xlabel('n'), 
title('ordinary least-squares solution, x(n), \lambda = 0'),
axis([0, 200, -5, 11], 'square'), text(120,9,'percent error = 192.78')


subplot(1,2,2)
plot(n,s,'b:',n,x1L2,'r'), xlabel('n'), 
title('L2 - regularized solution, x(n), \lambda = 0.01'),
axis([0, 200, -5, 11], 'square'), text(120,9,'percent error = 94.95'), 
legend('s(n)','x(n)','Location','southeast');
axis square

%% 1.c. L2-Regularized Solution, lambda = 0.1
lambda = 0.1;
x1L2 = (lambda*speye(L) + H'*H)\H'*y;

PerrorXL2 = 100*norm(x1L2 - s)/norm(s);

% Plots

figure(3)
plot(n,s,'b:',n,x1L2,'r'), xlabel('n'), 
title('L2 - regularized solution, x(n), \lambda = 0.1'),
axis([0, 200, -5, 11], 'square'), text(120,9,'percent error = 91.72'), 
legend('s(n)','x(n)','Location','southeast');


%% 1.d. L0- and L1-regularized Solutions Using CVX and IRLS
% L1 case - CVX solution
lambda = 0.1; p = 1; q = 1, eps = 1e-5, K = 100;

cvx_begin                                 
   variable x1CVXL1(L)
   minimize(sum_square(H*x1CVXL1-y) + lambda*norm(x1CVXL1,1))
cvx_end

Perr1CVXL1 = 100*norm(x1CVXL1-s)/norm(s);             % reconstruction error


% L1 case - IRLS solution

% Initialize IRLS algorithm with L2-regularized solution
x1L2 = (lambda*speye(L) + (H.')*H)\(H.'*y); 

x1IRLSL1 = x1L2;

for k = 1:K
WIRLS = diag(1./(abs(x1IRLSL1).^q + eps));
xtemp = (lambda*WIRLS + H.'*H)\(H.'*y);
Perr1IRLSL1(k) = 100*norm(xtemp - x1IRLSL1)/norm(xtemp);
x1IRLSL1 = xtemp;
end

PercentError1IRLSL1 = 100*norm(x1IRLSL1-s)/norm(s);

% L0 Case - CVX Solution
lambda = 0.1; p = 0; q = 2, eps = 1e-5, K = 100;

cvx_begin                                 
   variable x1CVXL0(L)
   minimize(sum_square(H*x1CVXL0-y) + lambda*norm(x1CVXL0,1))
cvx_end

Perr1CVXL0 = 100*norm(x1CVXL0-s)/norm(s);             % reconstruction error

% L0 Case - IRLS Solution

% Initialize IRLS algorithm with L2-regularized solution
x1IRLSL0 = x1L2;

for k = 1:K
WIRLS = diag(1./(abs(x1IRLSL0).^q + eps));
xtemp = (lambda*WIRLS + H.'*H)\(H.'*y);
Perr1IRLSL0(k) = 100*norm(xtemp - x1IRLSL0)/norm(xtemp);
x1IRLSL0 = xtemp;
end

PercentError1IRLSL0 = 100*norm(x1IRLSL0-s)/norm(s);

% Plots

figure(4)
subplot(2,2,1)
plot(n,x1CVXL1,'r'), title('L1 - CVX solution, x(n), \lambda = 0.1'),
xlabel('n'), text(120, 10, 'percent error = 130.1'),
axis([0, 200, -5, 11], 'square')

subplot(2,2,3)
plot(Ln,x1IRLSL1,'r'),xlabel('n'), 
title('L1 - IRLS solution, x(n), \lambda = 0.1'),
axis([0, 200, -5, 11], 'square'), text(120, 10,'percent error = 128.5'), 
axis square


subplot(2,2,2)
plot(n,x1CVXL0,'r'), title('L0 - CVX solution, x(n), \lambda = 0.1'),
xlabel('n'), text(120, 10, 'percent error = 130.1'),
axis([0, 200, -5, 11], 'square')


subplot(2,2,4)
plot(Ln,x1IRLSL0,'r'),xlabel('n'), 
title('L0 - IRLS solution, x(n), \lambda = 0.1'),
axis([0, 200, -5, 11], 'square'), text(120, 10,'percent error = 128.5'), 

figure(5)
subplot(1,2,1)
k = 1:K;
plot(k,Perr1IRLSL1,'r',k,Perr1IRLSL1,'b.'), title('L1 - IRLS iteration error, P(k)'),
xlabel('iterations, k'), ylabel('percent'), axis square, grid on

subplot(1,2,2)
plot(k,Perr1IRLSL0,'r',k,Perr1IRLSL0,'b.'), title('L0 - IRLS iteration error, P(k)'),
xlabel('iterations, k'), ylabel('percent'), axis square, grid on

%% 2. Sparse Signal Recovery
load Asy

%% 2.a. Input Signal and Noisy Observation
S = length(s);
n = 1:S;

figure(6)
subplot(1,2,1)
plot(n,s),title('sparse input, s(n)'), xlabel('n'), 
axis([0, 2000, -1.1, 1.1], 'square');
%
Y = length(y);
yPad = [y; zeros(S-Y, 1)];
subplot(1,2,2)
plot(n,yPad),title('observations, y(n)'), xlabel('n'), 
axis([0, 2000, -20, 20], 'square');
%% 2.b.Minimum Norm and L2-regularized Solutions

% L2-Minimum Norm Solution
x2MinNorm = pinv(A)*y;

PercentError1IRLSL0 = 100*norm(x2MinNorm-s)/norm(s);

% L2-regularized solution, lambda = 0.1
lambda = 0.1;
x2L2 = (lambda*speye(S) + A'*A)\(A'*y);

PerrorX2L2 = 100*norm(x2L2 - s)/norm(s);

% Plots

figure(7)
subplot(1,2,1)
plot(n,x2MinNorm,'r'),title('L2, minimum-norm solution, x(n)'), xlabel('n'), 
axis([0, 2000, -1.1, 1.1], 'square'), text(1400, 0.9,'percent error = 71.21%')


subplot(1,2,2)
plot(n,x2L2,'r'), xlabel('n'), 
title('L2 - regularized solution, x(n), \lambda = 0.1'),
axis([0, 2000, -1.1, 1.1], 'square'), text(1200, 0.9,'percent error = 71.21%'), 
axis square

%% 2.c. L0- and L1-regularized Solutions Using CVX and IRLS

% L1 case - CVX solution
lambda = 0.1; p = 1; q = 1, eps = 1e-6, K = 20;

cvx_begin                                 
   variable x2CVXL1(S)
   minimize(sum_square(A*x2CVXL1-y) + lambda*norm(x2CVXL1,1))
cvx_end

Perr2CVXL1 = 100*norm(x2CVXL1-s)/norm(s);             % reconstruction error

% L1 case - IRLS solution

% Initialize IRLS algorithm with L2-regularized solution

x2IRLSL1 = x2L2;

for k = 1:K
WIRLS = diag(1./(abs(x2IRLSL1).^q + eps));
xtemp = (lambda*WIRLS + A.'*A)\(A.'*y);
Perr2IRLSL1(k) = 100*norm(xtemp - x2IRLSL1)/norm(xtemp);
x2IRLSL1 = xtemp;
end

PercentError2IRLSL1 = 100*norm(x2IRLSL1-s)/norm(s);

% L0 Case - IRLS Solution
p = 0; q = 2;
% Initialize IRLS algorithm with L2-regularized solution
x2IRLSL0 = x2L2;

for k = 1:K
WIRLS = diag(1./(abs(x2IRLSL0).^q + eps));
xtemp = (lambda*WIRLS + A.'*A)\(A.'*y);
Perr2IRLSL0(k) = 100*norm(xtemp - x2IRLSL0)/norm(xtemp);
x2IRLSL0 = xtemp;
end

PercentError2IRLSL0 = 100*norm(x2IRLSL0-s)/norm(s);

% Plots

figure(8)
subplot(1,2,1)
plot(n,x2CVXL1,'r'), title('L1 - CVX solution, x(n), \lambda = 0.1'),
xlabel('n'), text(1200, 0.9, 'percent error = 2.54%'),
axis([0, 2000, -1.1, 1.1], 'square')

subplot(1,2,2)
plot(n,x2IRLSL1,'r'),xlabel('n'), 
title('L1 - IRLS solution, x(n), \lambda = 0.1'),
axis([0, 2000, -1.1, 1.1], 'square'), text(1200, 0.9,'percent error = 2.54%'), 

figure(9)
plot(n,x2IRLSL0,'r'),xlabel('n'), 
title('L0 - IRLS solution, x(n), \lambda = 0.1'),
axis([0, 2000, -1.1, 1.1], 'square'), text(1200, 0.9,'percent error = 0.63%'), 

figure(10)
subplot(1,2,1)
k = 1:K;
plot(k,Perr2IRLSL1,'r',k,Perr2IRLSL1,'b.'), title('L1 - IRLS iteration error, P(k)'),
xlabel('iterations, k'), ylabel('percent'), axis([0, 20, 0, 90], 'square'), grid on

subplot(1,2,2)
plot(k,Perr2IRLSL0,'r',k,Perr2IRLSL0,'b.'), title('L0 - IRLS iteration error, P(k)'),
xlabel('iterations, k'), ylabel('percent'), axis([0, 20, 0, 90], 'square'), grid on
