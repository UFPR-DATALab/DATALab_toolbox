function Xpar=paretopred_DTLab(Xpred,mX,sX)
% Xpar=paretopred_DTLab(Xpred,mX,sX);
% Does Pareto scaling on prediction X matrix based on calibration X matrix
%
% USAGE :
% Xas=parpred_DTLab(Xpred,Xcal);
%
% INPUT :
% Xpred : Prediction data matrix
% Xcal : Calibration data matrix
%
% OUTPUT :
% Xpar : Prediction X matrix transformed by Pareto scaling

[n,~]=size(Xpred);

% Mean center

Xmc = Xpred - ones(n,1) * mX;

% Pareto
Xpar= Xmc ./ (ones(n,1) * (sqrt(sX)));
