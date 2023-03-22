function Xmc=meancent_DTLab(X)
% Does Mean Center on X matrix
%
% USAGE :
% Xmc=meancent_DTLab(X);
%
% INPUT :
% X : Data matrix
%
% OUTPUT :
% Xmc : X matrix transformed by Mean Center

[n,~]=size(X);

% Mean center
mX = mean(X);
Xmc = X - ones(n,1) * mX;
