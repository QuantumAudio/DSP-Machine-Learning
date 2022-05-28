% p7ex.m  -  project 7 example 
% --------------------------------------
% DSPD - F 2017 - S. J. Orfanidis
% DSPD - S 2019 
%
% nicor.xls with PMA, and L1 and L2 W-H
% --------------------------------------

clear

Y = xlsread('nicor.xls');    % NICOR gas price data

y = Y(:,4);                  % closing prices
N = length(y);
t = 0:N-1; 
t7 = t/7;                    % t7 in weeks

% PMA/linear regression indicator
% -------------------------------

M = 19;                      % filter order
a = pma(y,M+1,0);            % filter with h0(n)  

% h = @(M) 2*(2*M+1-3*(0:M))/(M+1)/(M+2);   % equivalent calculation
% a = filter(h(M),1,y);
% a(1) = y(1);
% for K = 2:M,
%    hk = h(K-1);
%    yk = filter(hk,1,y(1:K));
%    a(K) = yk(end);
% end


% Whittaker-Henderson / L2
% ------------------------

s = 2; la = 20;         
N = length(y);          

Ds = diff(eye(N),s);
ts = t(s:end-1); ts7 = ts/7;

xw = (eye(N) + la*Ds'*Ds) \ y;

% Whittaker-Henderson / L1 with CVX
% ---------------------------------

la1 = 10;

cvx_quiet(true);
cvx_begin
    variable x1(N)
    minimize( sum_square(y-x1) + la1 * norm(Ds*x1,1) )
cvx_end

% Whittaker-Henderson / L0 with IRLS
% ----------------------------------

la0 = 1/2;

p = 0; q = 2 - p; 
ep = 1e-8; 
I = speye(N); 
K = 40; 

xold = (I + la0*Ds'*Ds) \ y;

for k=1:K,
   W = diag(1./(abs(Ds*xold).^q + ep));
   xk = (I + la0*Ds'*W*Ds) \ y;
   P(k) = 100 * norm(xk-xold)/norm(xold);
   xold = xk;
end

% --------------------------------------------

figure;
plot(t7,a,'g-', 'linewidth',2); hold on;
plot(t7,xw,'r-', 'linewidth',2);
ohlc(t7,Y);                          % uses O,H,L,C columns of Y
yaxis(41,51, 41:2:51)
xaxis(0,18,0:2:18); grid
title(['PMA,  {\itM} = ',num2str(M), '  /  L_2-WH,  {\its} = ', num2str(s), ', \lambda = ', num2str(la)])
xlabel('{\itt}  (weeks)');
legend(' pma', ' WH', ' data', 'location','se');

figure;
plot(t7,a,'g-', 'linewidth',2); hold on;
plot(t7,x1,'r-', 'linewidth',2);
ohlc(t7,Y); 
yaxis(41,51, 41:2:51)
xaxis(0,18,0:2:18); grid
title(['PMA,  {\itM} = ',num2str(M), '  /  L_1-WH,  {\its} = ', num2str(s), ', \lambda = ', num2str(la1)])
xlabel('{\itt}  (weeks)');
legend(' pma', ' WH', ' data', 'location','se');

figure;
plot(t7,a,'g-', 'linewidth',2); hold on;
plot(t7,xk,'r-', 'linewidth',2);
ohlc(t7,Y); 
yaxis(41,51, 41:2:51)
xaxis(0,18,0:2:18); grid
title(['PMA,  {\itM} = ',num2str(M), '  /  L_0-WH,  {\its} = ', num2str(s), ', \lambda = ', num2str(la0)])
xlabel('{\itt}  (weeks)');
legend(' pma', ' WH', ' data', 'location','se');


% L1 and L2 differenced signals, L1 is sparse
% -------------------------------------------

figure; plot(ts7,Ds*x1,'b-', ts7,Ds*xw, 'r--');
xaxis(0,18,0:2:18);
yaxis(-0.5,0.5, -0.5:0.25:0.5)
title(['\nabla^{\its}{\itx}({\itn}),  {\its} = ',num2str(s)])
xlabel('{\itt}  (weeks)');grid;
legend(' L_1', ' L_2', 'location','sw')

figure; plot(ts7,Ds*x1,'b-', ts7,Ds*xk, 'r--');
xaxis(0,18,0:2:18);
yaxis(-0.5,0.5, -0.5:0.25:0.5)
title(['\nabla^{\its}{\itx}({\itn}),  {\its} = ',num2str(s)])
xlabel('{\itt}  (weeks)');grid;
legend(' L_1', ' L_0', 'location','sw')

k = 1:K; 
figure; plot(k,P,'r-',k,P,'b.', 'markersize',16)
title('{\itL}_0,  IRLS iteration error');
ylabel('percent error'); xlabel('interation,  \itk');
yaxis(-0.1,0.4,0:0.1:0.4); grid 




















