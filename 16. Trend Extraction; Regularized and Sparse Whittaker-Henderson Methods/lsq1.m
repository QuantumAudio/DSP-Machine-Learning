% lsq1.m - least squares example

clear

A = [1 2 1 1; 1 1 -1 2]', x0 = [2 1]', b0 = A*x0

% A =
%      1     1
%      2     1
%      1    -1
%      1     2
% b0 =
%      3
%      5
%      1
%      4

b = [2.9, 5.2, 1.1, 3.9]'; 

x = A\b               

% x =
%     2.0879
%     0.9212

A'*A

%   7     4
%   4     7

sym((A'*A)^(-1))

%    [  7/33, -4/33]
%    [ -4/33,  7/33]

A'*b0
   
%    18.0000
%    15.0000

A'*b
   
%    18.3000
%    14.8000

x = (A'*A)^(-1) * (A'*b)

% x =
%     2.0879
%     0.9212

%% cvx -  http://cvxr.com/cvx/
% L1, L2 solutions

cvx_quiet(true)
cvx_begin
   variable x(2)
   minimize(norm(b-A*x,2))     % L2 norm
cvx_end

x

% x =
%     2.0879
%     0.9212


cvx_quiet(true)
cvx_begin
   variable x(2)
   minimize(norm(b-A*x,1))     % L1 norm
cvx_end

x

% x =
%     2.0333
%     0.9333

















