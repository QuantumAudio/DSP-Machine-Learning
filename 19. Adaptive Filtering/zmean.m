% zmean.m - zero mean of each column of a data matrix (or row vector)
%
% Usage: [Z,M] = zmean(Y)
%
% Y = N x p data matrix
% Z = N x p data matrix with zero-mean columns, or a zero-mean row
% M = N x p matrix whose columns are the means
%
% note: each column of Y is replaced by its zero-mean version,
%       Z = Y - ones(N,1) * m, where m = mean(Y)
%
%       Z = Y - M, or, Y = Z + M
%
%       if Y is a row vector, then Z is also a row vector
%
%       see also ZSTD
%
% S. J. Orfanidis - OSP Toolbox - 1999

function [Z,M] = zmean(Y)

if nargin==0, help zmean; return; end

M = ones(size(Y,1),1) * mean(Y);

Z = Y-M;


