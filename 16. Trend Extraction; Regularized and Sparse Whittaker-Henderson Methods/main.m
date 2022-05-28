%% 1a: Simulation Example 1a

% Load noisy version of flat-top signal
load yflat;
N = length(x);
n = 0:length(x) - 1;
% Plot exact input
close all
plot(n/2,x,'r');title('exact input x(n)');xlabel('n');grid on;
axis([0,300,0,5]);
% Plot noisy input
close all
plot(n/2,y,'b');title('noisy input y(n)');xlabel('n');grid on;
axis([0,300,0,5]);

%% 1b: Simulation Example 1b

% number of backward-difference iterations (derivatives)
s = 1; 
% regularization parameter
lambda = 5;
Ds = diff(speye(N),s);      % (N - s) x N, s-difference convolution matrix
x1b = (speye(N) + lambda*(Ds.')*Ds)\y; 
Dsx1b = Ds*x1b;

% Plot results
nDsx = 0:length(Dsx1b)-1;
close all
plot(n/2,x1b);title('L_2 criterion, s = 1, \lambda = 5');xlabel('n');
grid on;axis([0,300,0,5]);

close all
plot(nDsx/2,Dsx1b);title('L_2 case, \Delta^sx(n), s = 1, \lambda = 5');
xlabel('n');grid on;axis([0,300,-3, 3]);

%% 1c: Using the CVX Package
cvx_begin;
    variable x1cCVX(N);
    minimize(sum_square(x1cCVX - y) + lambda*norm(Ds*x1cCVX,1));
cvx_end;

Dsx1cCVX = Ds*x1cCVX;
% Plot Results
close all
plot(n/2,x1cCVX);title('L_1 criterion, CVX, s = 1, \lambda = 5');
xlabel('n');grid on;axis([0,300,0,5]);
%
close all
plot(nDsx/2,Dsx1cCVX);title('L_1 case, \Delta^sx(n), s = 1, \lambda = 5');
xlabel('n');grid on;axis([0,300,-3, 3]);

%% 1c: Using the IRLS L1 Method
K = 40;     % Number of iterations
% Inital x is defined as x1b in part 1b
x1cIRLSL1 = x1b;
p = 1; % L1 norm
q = 2-p; % To define denominator of weighted L2 norm
eps = 1e-10;

for k = 1:K;
W1cIRLS = diag(1./(abs(Ds*x1cIRLSL1).^q + eps));
xtemp = (speye(N) + lambda*(Ds.')*W1cIRLS*Ds)\y;
PL1(k) = 100*norm((xtemp - x1cIRLSL1), 2)/norm(xtemp, 2);
x1cIRLSL1 = xtemp;
end

Dsx1cIRLSL1 = Ds*x1cIRLSL1;
% Plot Results
close all
plot(n/2,x1cIRLSL1);title('L_1 criterion, IRLS, s = 1, \lambda = 5');
xlabel('n');grid on;axis([0,300,0,5]);
%
close all
plot(nDsx/2,Dsx1cIRLSL1);title('L_1 case, \Delta^sx(n), s = 1, \lambda = 5');
xlabel('n');grid on;axis([0,300,-3, 3]);
%
k = 1:K;
close all
plot(k, PL1,'r',k,PL1,'b.'),title('L_1, IRLS iteration error');xlabel('k');
ylabel('percent'),grid on;axis([0,40,-0.25,8])
%% 1c: Using the IRLS L0 Method
K = 40;     % Number of iterations
% Inital x is defined as x1b in part 1b
x1cIRLSL0 = x1b;
p = 0; % L0 norm
q = 2-p; % To define denominator of weighted L2 norm
eps = 1e-10;

for k = 1:K;
W1cIRLS = diag(1./(abs(Ds*x1cIRLSL0).^q + eps));
xtemp = (speye(N) + lambda*(Ds.')*W1cIRLS*Ds)\y;
PL0(k) = 100*norm((xtemp - x1cIRLSL0), 2)/norm(xtemp, 2);
x1cIRLSL0 = xtemp;
end

Dsx1cIRLSL0 = Ds*x1cIRLSL0;
% Plot Results
close all
plot(n/2,x1cIRLSL0);title('L_0 criterion, IRLS, s = 1, \lambda = 5');
xlabel('n');grid on;axis([0,300,0,5]);
%
close all
plot(nDsx/2,Dsx1cIRLSL0);title('L_0 case, \Delta^sx(n), s = 1, \lambda = 5');
xlabel('n');grid on;axis([0,300,-3, 3]);
%
k = 1:K;
close all
plot(k, PL0,'r',k,PL0,'b.'),title('L_0, IRLS iteration error');xlabel('k');
ylabel('percent'),grid on;axis([0,40,-0.25,8])

%% 1d: Comparing the EMA and SMA

% EMA
alpha = 0.95;
b = 1-alpha;
a = [1,-alpha];
x1cEMA = filter(b,a,y);
% Plot Results
close all
plot(n/2,x1cEMA),title('EMA output, \lambda = 0.95');xlabel('n');grid on;
axis([0,300,0,5]);

% SMA
L = round((1 + alpha)/(1 - alpha));
bSMA = (1/L)*ones(L,1);
% Plot Results
x1cSMA = filter(bSMA,1,y);
close all
plot(n/2,x1cSMA),title('SMA output, N = 39');xlabel('n');grid on;
axis([0,300,0,5]);

%% 2a: Simulation Example 2a

% Load noisy piecewise linear signal
load ylin;
N = length(x);
n = 0:length(x) - 1;
% Plot exact input
close all
plot(n/2,x,'r');title('exact input x(n)');xlabel('n');grid on;
axis([0,300,0,5]);
% Plot noisy input
close all
plot(n/2,y,'b');title('noisy input y(n)');xlabel('n');grid on;
axis([0,300,0,5]);

%% 2b: Simulation Example 2b

% number of backward-difference iterations (derivatives)
s = 2; 
% regularization parameter
lambda = 100;
Ds = diff(speye(N),s);      % (N - s) x N, s-difference convolution matrix
x2b = (speye(N) + lambda*(Ds.')*Ds)\y; 
Dsx2b = Ds*x2b;
% Plot Results
nDsx = 0:length(Dsx2b)-1;
close all
plot(n/2,x2b);title('L_2 criterion, s = 2, \lambda = 100');xlabel('n');
grid on;axis([0,300,0,5]);
%
close all
plot(nDsx/2,Dsx2b);title('L_2 case, \Delta^sx(n), s = 2, \lambda = 100');
xlabel('n');grid on;axis([0,300,-0.06, 0.06]);
%% 2c: Using the CVX Package
cvx_begin;
    variable x2cCVX(N);
    minimize(sum_square(x2cCVX - y) + lambda*norm(Ds*x2cCVX,1));
cvx_end;

Dsx2cCVX = Ds*x2cCVX;
% Plot Results
close all
plot(n/2,x2cCVX);title('L_1 criterion, CVX, s = 2, \lambda = 100');
xlabel('n');grid on;axis([0,300,0,5]);
%
close all
plot(nDsx/2,Dsx2cCVX);title('L_1 case, \Delta^sx(n), s = 2, \lambda = 100');
xlabel('n');grid on;axis([0,300,-0.06, 0.06]);

%% 2c: Using the IRLS L1 Method
K = 40;     % Number of iterations
% Inital x is defined as x2b in part 2b
x2cIRLSL1 = x2b;
p = 1; % L1 norm
q = 2-p; % To define denominator of weighted L2 norm
eps = 1e-10;

for k = 1:K
W2cIRLS = diag(1./(abs(Ds*x2cIRLSL1).^q + eps));
xtemp = (speye(N) + lambda*(Ds.')*W2cIRLS*Ds)\y;
PL1(k) = 100*norm((xtemp - x2cIRLSL1), 2)/norm(xtemp, 2);
x2cIRLSL1 = xtemp;
end

Dsx2cIRLSL1 = Ds*x2cIRLSL1;
% Plot Results
close all
plot(n/2,x2cIRLSL1);title('L_1 criterion, IRLS, s = 2, \lambda = 100');
xlabel('n');grid on;axis([0,300,0,5]);
%
close all
plot(nDsx/2,Dsx2cIRLSL1);title('L_1 case, \Delta^sx(n), s = 2, \lambda = 100');
xlabel('n');grid on;axis([0,300,-0.06, 0.06]);
%
k = 1:K;
close all
plot(k, PL1,'r',k,PL1,'b.'),title('L_1, IRLS iteration error');xlabel('k');
ylabel('percent'),grid on;axis([0,40,-0.25,4])
%% 2c: Using the IRLS L0 Method
lambda = 1;
K = 40;     % Number of iterations
% Inital x is defined as x1b in part 1b
x2cIRLSL0 = x2b;
p = 0; % L0 norm
q = 2-p; % To define denominator of weighted L2 norm
eps = 1e-10;

for k = 1:K
W2cIRLS = diag(1./(abs(Ds*x2cIRLSL0).^q + eps));
xtemp = (speye(N) + lambda*(Ds.')*W2cIRLS*Ds)\y;
PL0(k) = 100*norm((xtemp - x2cIRLSL0), 2)/norm(x2cIRLSL0, 2);
x2cIRLSL0 = xtemp;
end

Dsx2cIRLSL0 = Ds*x2cIRLSL0;
% Plot Results
close all
plot(n/2,x2cIRLSL0);title('L_0 criterion, IRLS, s = 2, \lambda = 1');
xlabel('n');grid on;axis([0,300,0,5]);
%
close all
plot(nDsx/2,Dsx2cIRLSL0);title('L_0 case, \Delta^sx(n), s = 2, \lambda = 1');
xlabel('n');grid on;axis([0,300,-0.06, 0.06]);
%
k = 1:K;
close all
plot(k, PL0,'r',k,PL0,'b.'),title('L_0, IRLS iteration error');xlabel('k');
ylabel('percent'),grid on;axis([0,40,-0.25,8])
%% 2d: Comparing the DEMA and SMA

% DEMA
[a,b,a1,a2,cinit] = dema(y,N); % Using output of first EMA as input to second EMA at time n = N - 1
% Plot Results
close all
%plot(n/2,x2cDEMA),title('DEMA output, \lambda = 0.95');xlabel('n');grid on;
axis([0,300,0,5]);

% SMA
L = round((1 + alpha)/(1 - alpha));
bSMA = (1/L)*ones(L,1);
x2cSMA = filter(bSMA,1,y);
% Plot Results
close all
plot(n/2,x2cSMA),title('SMA output, N = 39');xlabel('n');grid on;
axis([0,300,0,5]);

%% 3: Global Warming Trends

% Load global warming data
load tavenh2v.dat
% Extract dates
dates = tavenh2v(:,1); avg = tavenh2v(:,end);
N = length(avg);
n = 0:length(avg) - 1;
%% 3a: Using the Standard L2 Method
s = 2; lambda = 10000;
Ds = diff(speye(N),s);      % (N - s) x N, s-difference convolution matrix
x3 = (speye(N) + lambda*(Ds.')*Ds)\avg; 
Dsx3 = Ds*x3;
nDsx = linspace(dates(1),dates(end),length(Dsx3));
% Plot Results
close all
plot(dates,avg,'r',dates,x3,'b');title('L_2 case, s = 2, \lambda = 10,000');
xlabel('years');axis([1850,2000,-0.8,0.8]);legend('actual','trend');
ylabel('anomalies (C^o)');
%
close all
plot(nDsx,Dsx3);title('L_2 case, \Delta^sx(n), s = 2, \lambda = 10,000');
xlabel('n');grid on;axis([1850,2000,-0.03, 0.03]);

%% 3b: Using the CVX Package, L1 Method
lambda = 10;
cvx_begin;
    variable x3CVX(N);
    minimize(sum_square(x3CVX - avg) + lambda*norm(Ds*x3CVX,1));
cvx_end;

Dsx3CVX = Ds*x3CVX;
nDsx = linspace(dates(1),dates(end),length(Dsx3CVX));
% Plot Results
close all
plot(dates,avg,'r',dates,x3CVX,'b');title('L_1 case, CVX, s = 2, \lambda = 10');
xlabel('years');axis([1850,2000,-0.8,0.8]);legend('actual','trend');
ylabel('anomalies (C^o)');
%
close all
plot(nDsx,Dsx3CVX);title('L_1 case, CVX, \Delta^sx(n), s = 2, \lambda = 10');
xlabel('n');grid on;axis([1850,2000,-0.03, 0.03]);

%% 3c: Using the IRLS L1 Method
lambda = 10;
K = 40;     % Number of iterations
% Inital x is defined as x3 in the beginning of part 3
x3IRLSL1 = x3;
p = 1; % L1 norm
q = 2-p; % To define denominator of weighted L2 norm
eps = 1e-6;

for k = 1:K
W3IRLS = diag(1./(abs(Ds*x3IRLSL1).^q + eps));
xtemp = (speye(N) + lambda*(Ds.')*W3IRLS*Ds)\avg;
PL1(k) = 100*norm((xtemp - x3IRLSL1), 2)/norm(xtemp, 2);
x3IRLSL1 = xtemp;
end

Dsx3IRLSL1 = Ds*x3IRLSL1;
% Plot Results
close all
plot(dates,avg,'r',dates,x3IRLSL1,'b');title('L_1 case, IRLS, s = 2, \lambda = 10, \epsilon = 10^{-6}');
xlabel('years');axis([1850,2000,-0.8,0.8]);legend('actual','trend');
ylabel('anomalies (C^o)');
%
close all
plot(nDsx,Dsx3IRLSL1);title('L_1 case, IRLS, \Delta^sx(n), s = 2, \lambda = 10, \epsilon = 10^{-6}');
xlabel('n');grid on;axis([1850,2000,-0.03, 0.03]);
%
k = 1:K;
close all
plot(k, PL1,'r',k,PL1,'b.'),title('L_1 case, IRLS iteration error, \epsilon = 10^{-6}');xlabel('k');
ylabel('percent'),grid on;axis([0,40,-0.25, 20]);

%% 3: Using the IRLS L0 Method
lamda = 1;
K = 40;     % Number of iterations
% Inital x is defined as x3 in the beginning of part 3
x3IRLSL0 = x3;
p = 0; % L0 norm
q = 2-p; % To define denominator of weighted L2 norm
eps = 1e-3;

for k = 1:K
W3IRLS = diag(1./(abs(Ds*x3IRLSL0).^q + eps));
xtemp = (speye(N) + lambda*(Ds.')*W3IRLS*Ds)\avg;
PL0(k) = 100*norm((xtemp - x3IRLSL0), 2)/norm(xtemp, 2);
x3IRLSL0 = xtemp;
end

Dsx3IRLSL0 = Ds*x3IRLSL0;
% Plot Results
close all
plot(dates,avg,'r',dates,x3IRLSL0,'b');title('L_0 case, IRLS, s = 2, \lambda = 1, \epsilon = 10^{-3}');
xlabel('years');axis([1850,2000,-0.8,0.8]);legend('actual','trend');
ylabel('anomalies (C^o)');
%
close all
plot(nDsx,Dsx3IRLSL0);title('L_0 case, IRLS, \Delta^sx(n), s = 2, \lambda = 1, \epsilon = 10^{-3}');
xlabel('n');grid on;axis([1850,2000,-0.03, 0.03]);
%
k = 1:K;
close all
plot(k, PL0,'r',k,PL0,'b.'),title('L_0 case, IRLS iteration error, \epsilon = 10^{-3}');xlabel('k');
ylabel('percent'),grid on;axis([0,40,-0.25, 20]);