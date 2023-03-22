function Xas=autoscalpred_DTLab(Xpred,mX,sX)
% Does Autoscaling on prediction X matrix based on calibration X matrix
%
% USAGE :
% Xas=autoscalpred_DTLab(Xpred,Xcal);
%
% INPUT :
% Xpred : Prediction data matrix
% Xcal : Calibration data matrix
%
% OUTPUT :
% Xas : Prediction X matrix transformed by Autoscaling

[n,~]=size(Xpred);

% Mean center

Xmc = Xpred - ones(n,1) * mX;

% Autoscaling
Xas= Xmc ./ (ones(n,1) * sX);