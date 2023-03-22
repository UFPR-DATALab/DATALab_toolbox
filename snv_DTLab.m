function [Xsnv,mx,sx] = snv_DTLab(X)
% [Xsnv,mx,sx] = snv_DTLab(X)
% Does Standard Normal Variate Transformation (SNV) on X matrix
%
% USAGE :
% Xsnv=snv_DTLab(X);
%
% INPUT :
% X : Data matrix
%
% OUTPUT :
% Xsnv : X matrix transformed by snv
% mX : mean Data matrix
% sX : standard deviation Data matrix
%

[~,n]=size(X);
mx = mean(X,2);
sx = std(X,[],2);

Xsnv=(X-mx*ones(1,n))./(sx*ones(1,n));