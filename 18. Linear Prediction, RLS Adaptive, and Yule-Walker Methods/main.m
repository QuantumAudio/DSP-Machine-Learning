%% 1: AR(1) Model

Na = 500; Nb = 1500; N = 2000;

a1a = -0.75*ones(1, Na); % a1(n) for 0 <= n <= Na - 1
n1 = Na:1:Nb;             
a1b = -0.75*cos(pi*(n1 - Na)/(Nb - Na));  % a1(n) for Na <= n <= Nb
a1c = 0.75*ones(1,(N - Nb - 1));          % a1(n) for Nb + 1 <= n <= N - 1

a1 = [a1a, a1b, a1c]';

% White noise
rng('default');                     % Produces the same numbers each time
eps = normrnd(0, 1, [2000,1]);      % zero mean, unit variance

y(1) = 0;

% Synthesize signal
for n = 2:N
y(n, 1) = -a1(n)*y(n - 1, 1) + eps(n);
end

figure(1)
nPlot = 1:N;
plot(nPlot, y), title('AR(1) time-varying model, y(n)'), ylabel('y(n)'),
xlabel('time samples, n')

% Analysis Equations, la = 0.98
la = 0.98; delta = 1e-4; R0(1) = delta; R1(1) = 0; a1hat(1) = 0;

for n = 2:N
    R0(n) = la*R0(n - 1) + y(n - 1)^2;
    R1(n) = la*R1(n - 1) + y(n)*y(n - 1);
    a1hat(n) = -R1(n)/R0(n);
end

figure (2)
plot(nPlot,a1hat,'b', nPlot, a1, 'r--'), title('AR(1) parameter tracking, \lambda = 0.980')
grid on, legend('a_1(n) adaptive', 'a_1(n) theoretical'), xlabel('time samples, n')

% Analysis Equations, la = 0.995
la = 0.995; delta = 1e-4; R0(1) = delta; R1(1) = 0; a1hat(1) = 0;

for n = 2:N
    R0(n) = la*R0(n - 1) + y(n - 1)^2;
    R1(n) = la*R1(n - 1) + y(n)*y(n - 1);
    a1hat(n) = -R1(n)/R0(n);
end


figure (3)
plot(nPlot,a1hat,'b', nPlot, a1, 'r--'), title('AR(1) parameter tracking, \lambda = 0.995')
grid on, legend('a_1(n) adaptive', 'a_1(n) theoretical'), xlabel('time samples, n')

%% 2a: AR(2) Model

% Na = 500; Nb = 1500; N = 2000;

a1a = -1.3*ones(1, Na); % a1(n) for 0 <= n <= Na - 1
n1 = Na:1:Nb;             
a1b = 1.3*(n1 - Nb)/(Nb - Na);      % a1(n) for Na <= n <= Nb
a1c = zeros(1,(N - Nb - 1));    % a1(n) for Nb + 1 <= n <= N - 1

a1 = [a1a, a1b, a1c]';

a2a = 0.4*ones(1, Na); % a1(n) for 0 <= n <= Na - 1
n1 = Na:1:Nb;             
a2b = 0.65 - 0.25*cos(pi*(n1 - Na)/(Nb - Na));      % a1(n) for Na <= n <= Nb
a2c = 0.9*ones(1,(N - Nb - 1));    % a1(n) for Nb + 1 <= n <= N - 1

a2 = [a2a, a2b, a2c]';

y(1) = 0;
y(2) = 0;

% Synthesize signal
for n = 3:N
y(n, 1) = -a1(n)*y(n - 1, 1) -a2(n)*y(n - 2, 1) + eps(n);
end

figure(4)
plot(nPlot, y), title('AR(2) time-varying model, y(n)'), ylabel('y(n)'),
xlabel('time samples, n')

%% 2.b: Computing Adaptive Coefficients

% Analysis Equations,la = 0.98
la = 0.98; delta = 1e-3; % using delta = 1e-4 gave NaN for at least one value

R11(1) = delta; R12(1) = delta; R22(1) = delta; R01(1) = delta;
R02(1) = delta; a1hat(1) = 0; a2hat(1) = 0;

R11(2) = delta; R12(2) = delta; R22(2) = delta; R01(2) = delta;
R02(2) = delta; a1hat(2) = 0; a2hat(2) = 0;

for n = 3:N

R11(n) = la*R11(n-1) + y(n-1)*y(n-1);
R12(n) = la*R12(n-1) + y(n-1)*y(n-2);
R22(n) = la*R22(n-1) + y(n-2)*y(n-2);
R01(n) = la*R01(n-1) + y(n)*y(n-1);
R02(n) = la*R02(n-1) + y(n)*y(n-2);

a12 = -[R11(n), R12(n); R12(n), R22(n)]\[R01(n); R02(n)];

a1hat(n) = a12(1); a2hat(n) = a12(2);
end

figure(5)
plot(nPlot, a1hat, nPlot, a1, 'r--'), grid on, xlabel('time samples, n')
title('AR(2) parameter tracking, \lambda = 0.980'),
legend('a_1(n) adaptive', 'a_1(n) theoretical')

figure(6)
plot(nPlot, a2hat, nPlot, a2, 'r--'), grid on, xlabel('time samples, n')
title('AR(2) parameter tracking, \lambda = 0.980'),
legend('a_2(n) adaptive', 'a_2(n) theoretical'), axis([0,2000,0,1.2]);

% Analysis Equations,la = 0.995
la = 0.995; delta = 1e-3; % using delta = 1e-4 gave NaN for at least one value

R11(1) = delta; R12(1) = delta; R22(1) = delta; R01(1) = delta;
R02(1) = delta; a1hat(1) = 0; a2hat(1) = 0;

R11(2) = delta; R12(2) = delta; R22(2) = delta; R01(2) = delta;
R02(2) = delta; a1hat(2) = 0; a2hat(2) = 0;

for n = 3:N

R11(n) = la*R11(n-1) + y(n-1)*y(n-1);
R12(n) = la*R12(n-1) + y(n-1)*y(n-2);
R22(n) = la*R22(n-1) + y(n-2)*y(n-2);
R01(n) = la*R01(n-1) + y(n)*y(n-1);
R02(n) = la*R02(n-1) + y(n)*y(n-2);

a12 = -[R11(n), R12(n); R12(n), R22(n)]\[R01(n); R02(n)];

a1hat(n) = a12(1); a2hat(n) = a12(2);
end

figure(7)
plot(nPlot, a1hat, nPlot, a1, 'r--'), grid on, xlabel('time samples, n')
title('AR(2) parameter tracking, \lambda = 0.995'),
legend('a_1(n) adaptive', 'a_1(n) theoretical'), axis([0,2000,-1.7,0.7]);

figure(8)
plot(nPlot, a2hat, nPlot, a2, 'r--'), grid on, xlabel('time samples, n')
title('AR(2) parameter tracking, \lambda = 0.995'),
legend('a_2(n) adaptive', 'a_2(n) theoretical'), axis([0,2000,0,1.2]);

%% 2.c AR(2) Analysis with AR(1) Synthesis

a1a = -0.75*ones(1, Na); % a1(n) for 0 <= n <= Na - 1
n1 = Na:1:Nb;             
a1b = -0.75*cos(pi*(n1 - Na)/(Nb - Na));  % a1(n) for Na <= n <= Nb
a1c = 0.75*ones(1,(N - Nb - 1));          % a1(n) for Nb + 1 <= n <= N - 1

a1 = [a1a, a1b, a1c]';

a2 = zeros(2000, 1);

y(1) = 0;
y(2) = 0;

% AR(1) Synthesis
for n = 3:N
y(n, 1) = -a1(n)*y(n - 1, 1) + eps(n);
end

% AR(2) Analysis
la = 0.98; delta = 1e-3; % using delta = 1e-4 gave NaN for at least one value

R11(1) = delta; R12(1) = delta; R22(1) = delta; R01(1) = delta;
R02(1) = delta; a1hat(1) = 0; a2hat(1) = 0;

R11(2) = delta; R12(2) = delta; R22(2) = delta; R01(2) = delta;
R02(2) = delta; a1hat(2) = 0; a2hat(2) = 0;

for n = 3:N

R11(n) = la*R11(n-1) + y(n-1)*y(n-1);
R12(n) = la*R12(n-1) + y(n-1)*y(n-2);
R22(n) = la*R22(n-1) + y(n-2)*y(n-2);
R01(n) = la*R01(n-1) + y(n)*y(n-1);
R02(n) = la*R02(n-1) + y(n)*y(n-2);

a12 = -[R11(n), R12(n); R12(n), R22(n)]\[R01(n); R02(n)];

a1hat(n) = a12(1); a2hat(n) = a12(2);
end

figure(9)
plot(nPlot, a1hat, nPlot, a1, 'r--'), grid on, xlabel('time samples, n')
title('AR(1) with order-2, \lambda = 0.980'),
legend('a_1(n) adaptive', 'a_1(n) theoretical'), axis([0, 2000, -1, 1]);

figure(10)
plot(nPlot, a2hat, nPlot, a2, 'r--'), grid on, xlabel('time samples, n')
title('AR(1) with order-2, \lambda = 0.980'),
legend('a_2(n) adaptive', 'a_2(n) theoretical'), axis([0, 2000, -1, 2]);

%% 3: Sunspot Data
%% 3a: Yule-Walker Method, order 2
Y = load('sunspots.dat', '-ascii');
i = find(Y(:,1) == 1809);       % Find index of year 1809
y = Y(i:end, 2);                % Extract all years since 1809
N = length(y);                  % N = 200
m = mean(y); y = y-m;           % Make data zero-mean

% compute lag-0 value
S = 0;    % temp variable
for n = 1:N
    temp = y(n)*y(n);
    S = S + temp;
end

R0hat = S/N;        

% compute lag-1 value
S = 0;    % temp variable
for n = 1:N-1
    temp = y(n+1)*y(n);
    S = S + temp;
end

R1hat = S/N;

% compute lag-2 value
S = 0;    % temp variable
for n = 1:N-2
    temp = y(n+2)*y(n);
    S = S + temp;
end

R2hat = S/N;

a12 = -[R0hat, R1hat; R1hat, R0hat]\[R1hat; R2hat];

a1hatEst = a12(1); a2hatEst = a12(2);
%a1hatEst = -1.3777; a2hatEst = 0.6822;

%% 3.b: RLS Method, order-2
clear a1hat; clear a2hat;
la = 0.999; delta = 1e-4; % using delta = 1e-4 gave NaN for at least one value

R11(1) = delta; R12(1) = delta; R22(1) = delta; R01(1) = delta;
R02(1) = delta; a1hat(1) = 0; a2hat(1) = 0;

R11(2) = delta; R12(2) = delta; R22(2) = delta; R01(2) = delta;
R02(2) = delta; a1hat(2) = 0; a2hat(2) = 0;

for n = 3:N

R11(n) = la*R11(n-1) + y(n-1)*y(n-1);
R12(n) = la*R12(n-1) + y(n-1)*y(n-2);
R22(n) = la*R22(n-1) + y(n-2)*y(n-2);
R01(n) = la*R01(n-1) + y(n)*y(n-1);
R02(n) = la*R02(n-1) + y(n)*y(n-2);

a12 = -[R11(n), R12(n); R12(n), R22(n)]\[R01(n); R02(n)];

a1hat(n) = a12(1); a2hat(n) = a12(2);
end

figure(11)
nPlot = 1:N;
plot(nPlot, a1hat, nPlot, a1hatEst*ones(1,N),'r--')
title('a_1 coefficient, \lambda = 0.999'), xlabel('years'),
legend('adaptive', 'Yule-Walker'), 

figure(12)
plot(nPlot, a2hat, nPlot, a2hatEst*ones(1,N),'r--')
title('a_2 coefficient, \lambda = 0.999'), xlabel('years'),
legend('adaptive', 'Yule-Walker'), 

%% 3.c: AR(2) Synthesis

% RLS Method
yAdapt(1) = 0;
yAdapt(2) = 0;

for n = 3:N
    yAdapt(n) = -a1hat(n)*y(n-1) -a2hat(n)*y(n-2);
end

yAdapt = yAdapt + m;    % add back mean value

% Yule-Walker Method
yYW(1) = 0;
yYW(2) = 0;

for n = 3:N
    yYW(n) = -a1hatEst*y(n-1) -a2hatEst*y(n-2);
end

yYW = yYW + m;          % add back mean value

y = y + m;              % add back mean value

figure(13)
plot(nPlot, y, 'b:', nPlot, yAdapt, 'r', nPlot, yYW, 'g'), 
title('Sunspot Numbers 1809 - 2008'),
xlabel('years'), legend('data', ' adaptive prediction', 'Yule-Walker prediction')

%% 3.d: Comparing Spectra

% make data zero-mean again
y = y - m;

% Periodogram
n = 0:N-1;
% calculate Hamming window 
w = (0.54 - 0.46*cos(2*pi*n/N))';

p = 2:0.01:20;
omega = 2*pi./p;


[Sper, SperW] = periodogram(y,w,omega);
Sper = Sper/max(Sper);
Sper = [zeros(1, 2000 - length(Sper)), Sper]; % pad zeros to fill in 0 <= omega < 2

% Calculate AR(2) periodogram spectrum

variance = R0hat + a1hatEst*R1hat + a2hatEst*R2hat;

Sar = variance./abs(1 + a1hatEst.*exp(-j*omega) + a2hatEst.*exp(-2*j*omega)).^2;

% add values from 0 to 2 to graphs by padding data
Sar = [zeros(1, 199), Sar];
% normalize
Sar = Sar/max(Sar);


%% 3.e.I: RLS method, order-9, la = 0.999

Y = load('sunspots.dat', '-ascii');
i = find(Y(:,1) == 1809);       % Find index of year 1809
y = Y(i:end, 2);                % Extract all years since 1809
N = length(y);                  % N = 200
m = mean(y); y = y-m;           % Make data zero-mean

la = 0.999; delta = 1e-4; P = (1/delta)*speye(9); ybuff(1:9,1) = 0; a(1:9,1) = 0;
N = length(y);

% RLS Method
for n = 0:N-1
    k = (1/la)*P*ybuff;
    mu = 1/(1 + k'*ybuff);
    e = y(n+1) + a'*ybuff;
    a = a - mu*e*k;
    P = (1/la)*P - mu*(k*k');
    ybuff = [y(n+1); ybuff(1:(end - 1))];
    yAdapt(n+1) = -a'*ybuff;
end

% add back the mean

yAdapt = yAdapt + m;

%% 3.e.II: Yule-Walker method, order-9
Y = load('sunspots.dat', '-ascii');
i = find(Y(:,1) == 1809);       % Find index of year 1809
y = Y(i:end, 2);                % Extract all years since 1809
N = length(y);                  % N = 200
m = mean(y); y = y-m;           % Make data zero-mean

Ryy = acf(y,y);

R1 = [Ryy(1) Ryy(2) Ryy(3) Ryy(4) Ryy(5) Ryy(6) Ryy(7) Ryy(8) Ryy(9);
      Ryy(2)  Ryy(1) Ryy(2) Ryy(3) Ryy(4) Ryy(5) Ryy(6) Ryy(7) Ryy(8);
      Ryy(3)  Ryy(2) Ryy(1) Ryy(2) Ryy(3) Ryy(4) Ryy(5) Ryy(6) Ryy(7);
      Ryy(4)  Ryy(3) Ryy(2) Ryy(1) Ryy(2) Ryy(3) Ryy(4) Ryy(5) Ryy(6);
      Ryy(5)  Ryy(4) Ryy(3) Ryy(2) Ryy(1) Ryy(2) Ryy(3) Ryy(4) Ryy(5);
      Ryy(6)  Ryy(5) Ryy(4) Ryy(3) Ryy(2) Ryy(1) Ryy(2) Ryy(3) Ryy(4);
      Ryy(7)  Ryy(6) Ryy(5) Ryy(4) Ryy(3) Ryy(2) Ryy(1) Ryy(2) Ryy(3);
      Ryy(8)  Ryy(7) Ryy(6) Ryy(5) Ryy(4) Ryy(3) Ryy(2) Ryy(1) Ryy(2);
      Ryy(9)  Ryy(8) Ryy(7) Ryy(6) Ryy(5) Ryy(4) Ryy(3) Ryy(2) Ryy(1)];
    
R2 = [Ryy(2); Ryy(3); Ryy(4); Ryy(5); Ryy(6); Ryy(7); Ryy(8); Ryy(9); Ryy(10)];

aEst = -R1\R2;

ybuff(1:9,1) = 0;

for n = 0:N-1
    yYW9(n+1) = -aEst'*ybuff;
    ybuff = [y(n+1); ybuff(1:(end-1))];
end

yYW9 = yYW9 + m ; % add back mean
y = y + m; % add back mean
%

nPlot = 1:length(y);
figure(14)
plot(nPlot, y, 'b:', nPlot, yAdapt, 'r', nPlot, yYW9, 'g'), 
title('Sunspot Numbers 1809 - 2008'),
xlabel('years'), legend('data', ' adaptive prediction', 'Yule-Walker prediction')
axis([0, 200, 0, 200]);
% Calculate AR(9) periodogram spectrum

RTemp = Ryy(1:10)';
aEstTemp = [1; aEst]; % add a coefficient of 1 to multiply by R0hat
variance9 = aEstTemp'*RTemp;

% p = 2:0.01:20;        % Defined previously
% omega = 2*pi./p;

Sar9 = variance9./abs(1 + aEst(1)*exp(-j*omega) + ...
    aEst(2)*exp(-2*j*omega) + aEst(3)*exp(-3*j*omega) + ...
    aEst(4)*exp(-4*j*omega) + aEst(5)*exp(-5*j*omega) + ...
    aEst(6)*exp(-6*j*omega) + aEst(7)*exp(-7*j*omega) + ...
    aEst(8)*exp(-8*j*omega) + aEst(9)*exp(-9*j*omega)).^2;

% add values from 0 to 2 to graphs by padding data

Sar9 = [zeros(1, 199), Sar9];
% normalize
Sar9 = Sar9/max(Sar9);


% Plot periodogram
nPlot = 1:length(Sper);
nPlot2 = 1:length(Sar);
nPlot3 = 1:length(Sar9);

figure(15)
plot(nPlot*0.01, Sper, 'r:',nPlot2*0.01, Sar, 'k:', nPlot3*0.01, Sar9, 'b'); 
legend('periodogram', 'YW, AR(2)','YW, AR(9)'), grid on,
title('Periodogram and Yule-Walker Spectra'), xlabel('period in years'),
axis([0, 20, 0, 1.2])

