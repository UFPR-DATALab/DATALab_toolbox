function Xmc=meancentpred_DTLab(Xpred,mX)
% Does Mean center on prediction X matrix based on calibration X matrix
%
% USAGE :
% Xmc=meancentpred_DTLab(Xpred,mX)
%
% INPUT :
% Xpred : Prediction data matrix
% mx : Mean of Calibration data matrix
%
% OUTPUT :
% Xmc : Prediction X matrix transformed by mean center

[n,~]=size(Xpred);

% Mean center

Xmc = Xpred - ones(n,1) * mX;
