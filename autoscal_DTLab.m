function Xas=autoscal_DTLab(X)
% Does Autoscaling on X matrix
%
% USAGE :
% Xas=autoscal_DTLab(X);
%
% INPUT :
% X : Data matrix
%
% OUTPUT :
% Xas : X matrix transformed by Autoscaling

[n,~]=size(X);

% Mean center
mX = mean(X);
Xmc = X - ones(n,1) * mX;

% Autoscaling
Xas= Xmc ./ (ones(n,1) * std(Xmc));