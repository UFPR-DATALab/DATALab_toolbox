function [Model_pred] = plsdapred_DTLab(Model_pred,Xval,varargin)
% [Model_pred] = plsdapred_DTLab(Model,Xval,options,Yval);
%
% INPUT :
% Model : Output from plsdacal_DTLab
% Xval : Prediction/Validation X
% options.plots = display plot (0 (none) / 1 (some - DEFAULT) / 2 (all))
% Yval : Validation Y (For Blind test)
%
%
% OUTPUT:
% Model_pred
% Model_pred.Ypred : Predicted responses
% Model_pred.DModXn_pred: Normalized predicted DModX
% Model_pred.T2_pred: Predicted Hotteling
% Model_pred.RMSEV  = RMSEV for blind test
% Model_pred.R2V  = R²V for blind test
%
% version 1.0
% DATALab - Data Science in Chemistry and Chemometrics Laboratory
% 

%% Check options
% if exist('options','var')
%     if isfield(options,'plots')
%         plots = options.plots;
%     else
%         plots = 1;
%     end
% else
%     plots = 1;
% end

%% Model data recovery
ppx = Model_pred.ppx;
ppy = Model_pred.ppy;
LV = Model_pred.LV;
b = Model_pred.beta;
Ppred = Model_pred.Xloading;
Yc = Model_pred.Yc;
for i=1:size(Model_pred.stats,2)
    W{i,1} = Model_pred.stats{1,i}.W;
    Tcrit(i,1) = Model_pred.stats{1,i}.Tcrit;
    Qcrit(i,1) = Model_pred.stats{1,i}.Qcrit;
    Q(:,i) = Model_pred.stats{1,i}.Q;
end    
XS = Model_pred.Xscores;
if exist('Model.Threshold','var')
    Thrs = Model_pred.Threshold;
else
    Thrs = 0.5;
end



%% Pre-treatment
if ppx == 1
    Xpp_v = meancentpred_DTLab(Xval,Model_pred.mX);
elseif ppx == 2
    Xpp_v = autoscalpred_DTLab(Xval,Model_pred.mX,Model_pred.sX);
elseif ppx == 3
    Xpp_v = paretopred_DTLab(Xval,Model_pred.mX,Model_pred.sX);
else
    Xpp_v = Xval;
end

[uu] = size(Yc,2);
%% Prediction
Ypred = [ones(size(Xpp_v,1),1),Xpp_v]*b;

if ppy == 1
    Yp = Ypred + mean(Yc);
elseif ppy == 2
    Yp = Ypred*std(Yc) + mean(Yc);
elseif ppy == 3
    Yp = Ypred*sqrt(std(Yc)) + mean(Yc);
else
    Yp = Ypred;
end

%% DModX
for i=1:size(Model_pred.stats,2)
    Tpred{i,1} = Xpp_v*W{i,1};
    Xrec{i,1} = Tpred{i,1}*Ppred{1,i}';
    ex{i,1} = Xpp_v-Xrec{i,1};
    t2_pred{i,1} = sum( bsxfun(@rdivide, abs(Tpred{i,1}).^2, var(XS{1,i},[],1)) , 2);
    
    err = ex{i,1};
    for j=1:size(ex,1)
        Qp(i,1) = err(j,:)*err(j,:)';
    end    
    Q_pred{i,1} = Qp;
end

%% Confusion Matrix
Model_pred.Ypred = Yp;

% [~,Ypr] = max(Yp,[],2);
Ypr = Yp;
for j=1:size(Thrs,1)
    chck = Ypr(:,j)<Thrs(j);
    Ypr(chck,j) = 0;
end

% Return Groups to Y
Cls = zeros(size(Yp,1),1);
[~,ind] = max(Yp,[],2);
sCls = find(sum(Ypr,2)~=1);

for j=1:1:uu
    Cls(Ypr(:,j)==1) = j;
end
Cls(sCls) = ind(sCls);
Cls(sum(Ypr,2)==0) = 0;

if size(varargin,2) > 1
    Yval = varargin{1,2};
    d_Yp = Cls - Yval;
    Clas_errP = sum(d_Yp~=0,1)./sum(Yval,1);
    p_conf = confusionmat(Yval,Cls);
    p_conf = p_conf';
end

Model_pred.Class = Cls;
Model_pred.Q_pred = Q_pred;
Model_pred.T2_pred = t2_pred;
if size(varargin,2) > 1
    Model_pred.ClassficationError  = Clas_errP;
    Model_pred.ConfusionM_p  = p_conf;
end

%% Figures of Merit
if size(varargin,2) > 1
    Accuracy_p = (trace(p_conf)/(sum(p_conf(:))))*100;
    MissClassRate_p = 100-Accuracy_p;
    Sen_p = diag(p_conf)./sum(p_conf,2);
    Spec_p = (repmat(sum(p_conf(:)),[size(p_conf),1])-sum(p_conf,2)-sum(p_conf,1)+diag(p_conf))./(repmat(sum(p_conf(:)),[size(p_conf),1]) - sum(p_conf,2));
    
    Model_pred.CM_stats.Accuracy_p               = Accuracy_p;
    Model_pred.CM_stats.MissClassificationRate_p = MissClassRate_p;
    Model_pred.CM_stats.Sensibility_p            = Sen_p;
    Model_pred.CM_stats.Specificity_p            = diag(Spec_p);
    
%     % ROC curve
%     u = unique(Yval);
%     [uu,~] = size(u);
%     f = cell(uu,1);
%     for j=1:1:uu
%         f{j} = find(Yval == u(j,1));
%     end
%     Grp = zeros(size(Yval,1),uu);
%     for j=1:1:uu
%         Grp(f{j},j) = 1;
%     end
%     
%     Xroc = zeros(size(Ypr,2)+1,size(u,1));
%     Yroc = zeros(size(Ypr,2)+1,size(u,1));
%     AUC = zeros(size(u,1),1);
%     for j=1:1:uu
%         [Xroc(:,j),Yroc(:,j),~,AUC(j)] = perfcurve(Grp(:,j),Ypr(:,j),1);
%     end
%     
%     Model_pred.ROC.FPrate = Xroc;
%     Model_pred.ROC.TPrate = Yroc;
%     Model_pred.ROC.AUC    = AUC;
end
