function [VIP,VIP_ord,Ind]=vip_DTLab(X,XS,W,YL)
% Does VIP to variable selection
%
% USAGE :
% [VIP,VIP_ord,Ind]=vip_DTLab(X,XS,W,YL);
%
% INPUT :
% X : Data matrix
% Y : Response matrix
% XS : scores of X matrix from PLS method
% W : weight matrix from PLS method (stats.W from plsregress)
% YL : loadings of Y matrix from PLS method
%
% OUTPUT :
% VIP : vector with variable projection importance
% VIP_ord : Sorted VIP from highest to lowest VIP
% Ind : Variable index of Sorted VIP
%Reference: A review of variable selection methods in Partial Least Squares Regression
% Tahir Mehmood et al, Chemometrics and Intelligent Laboratory Systems 118 (2012)62-69
if size(YL,1)>1
    for k=1:1:size(YL,1)
    yl = YL(k,:);
    SS(:,k)=diag((yl*yl')*XS'*XS);
    end
else
    SS=diag((YL*YL')*XS'*XS);
end

%initializing
[~,n]=size(X);

%calculate VIP;
nW = vecnorm(W);
% nW = norm(W);
w = W./(ones(n,1)*nW);
w = w.^2;
expl_var=SS'*w';  % explained variance by variable i
VIP=sqrt((n*expl_var)'/sum(SS,1)); %Matriz de covariancia entre T e Y ponderada pelos pesos de cada LV, dividida pela matriz de covariância entre T e Y


if nargin > 1
    [VIP_ord,Ind] = sort(VIP,'ascend');
end

