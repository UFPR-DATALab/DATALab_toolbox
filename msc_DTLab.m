function [Xmsc]= msc_DTLab(X, Xref)
% [Xmsc]= msc_LQQA(X, Xref);
% Does Multiplicative Scattering Correction on X matrix
%
% USAGE :
% [Xmsc]= msc_DTLab(X, Xref);
%
% INPUT :
% X : Data matrix
% Xref : Reference spectrum for correction (Default is mean of X)
%
% OUTPUT :
% Xmsc : X matrix corrected by Multiplicative Scattering

if(nargin==1)
   Xref=mean(X);
end

[~,p]=size(X);
Xm = [ones(1,p);Xref];
coef=pinv(Xm')*X';
Xmsc=((X'-ones(p,1)*coef(1,:))./(ones(p,1)*coef(2,:)))';


