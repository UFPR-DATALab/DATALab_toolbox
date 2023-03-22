function [Model] = plsdacal_DTLab(Xcal,Ycal,LV,options)
% [Model] = plsdacal_DTLab(Xcal,Ycal,LV,options);
%
% INPUT :
% Xcal : Calibration X
% Ycal : Calibration Y
% LV : Maximum latent variables to perform the model
% options.threshold : Define the threshold method (1 = 0.5 / 2 = Baysean)
% options : PLS options - for more information: help plscal_DTLab
%
% OUTPUT:
% Model
% Model.Unclassified_cal = Unclassified samples during calibration;
% Model.Unclassified_cv = Unclassified samples during cross-validation;
% Model.ConfusionM_cal = Confusion matrix of calibration;
% Model.ConfusionM_cv = Confusion matrix of cross-validation;
% Model.stats : With Accuracy, Specificity, Selectivity and Miss
% Classification Rate
%
% USAGE:
% plsregress - Statistics and Machine Learning Toolbox
% plscal_DTLab - DATALab Toolbox
%
% version 1.0
% DATALab - Data Science in Chemistry and Chemometrics Laboratory
% 
%% Ycal to Groups
u = unique(Ycal);
[uu,~] = size(u);
f = cell(uu,1);
for j=1:1:uu
    f{j} = find(Ycal == u(j,1));
end
Grp = zeros(size(Xcal,1),uu);
for j=1:1:uu
    Grp(f{j},j) = 1;
end

%% PLS
options.ppy = 0;
options.pls_type = 2;
[Model] = plscal_DTLab(Xcal,Grp,LV,options);

if exist('options','var')
    if isfield(options,'threshold')
        trshld_mthd = options.threshold;
    else
        trshld_mthd = 1;
    end
else
    trshld_mthd = 1;
end


if trshld_mthd == 1 %0.5 threshold
    trshld = repmat(0.5,uu,1);
    yfit = round(Model.Yfit); %0.5 threshold
    ycv = round(Model.Ycv);
    Y = Model.Yc;
    Model.Threshold = trshld;
elseif trshld_mthd == 2 % Bayes approach
    yfit = Model.Yfit;
    ycv = Model.Ycv;
    Y = Model.Yc;
    options.plots = 0;
    trshld = zeros(uu,1);
    for j=1:uu
        [trshld(j,1), ~] = plsda_bayseatnrshld_DTLab(Y(:,j),yfit(:,j),options);
    end
    Model.Threshold = trshld;
end

% Return Groups to Y
Yc = zeros(size(yfit,1),1);
Ycv = zeros(size(ycv,1),1);
Ycal_new = zeros(size(ycv,1),1);
for j=1:1:uu
    Yc(yfit(:,j)>trshld(j)) = j;
    Ycv(ycv(:,j)>trshld(j)) = j;
    Ycal_new(Y(:,j)>trshld(j)) = j;
end

%% Confusion matrix
if size(unique(Yc),1)~=size(unique(Ycal_new),1)
   unclsGrp = setdiff(unique(Yc,'legacy'),unique(Ycal_new,'legacy'));
   cal_conf = confusionmat(Ycal_new,Yc);
   unclsC = cal_conf(end-(length(unclsGrp)-1):end,:);
   cal_conf(end-(length(unclsGrp)-1):end,:) = [];
   Model.Unclassified_cal = unclsC;
else
    cal_conf = confusionmat(Ycal_new,Yc);
    cal_conf = cal_conf';
end

if size(unique(Ycv),1)~=size(unique(Ycal_new),1)
   unclsGrp = setdiff(unique(Ycv),unique(Ycal_new));
   cv_conf = confusionmat(Ycal_new,Ycv);
   unclsCV = cv_conf(end-(length(unclsGrp)-1):end,:);
   cv_conf(end-(length(unclsGrp)-1):end,:) = [];
   Model.Unclassified_cv = unclsCV;
else
    cv_conf = confusionmat(Ycal_new,Ycv);
    cv_conf = cv_conf';
end


Model.ConfusionM_cal = cal_conf;
Model.ConfusionM_cv  = cv_conf;


%% Figures of Merit
Accuracy_cal = (trace(cal_conf)/(sum(cal_conf(:))))*100;
Accuracy_cv = (trace(cv_conf)/(sum(cv_conf(:))))*100;
MissClassRate_cal = 100-Accuracy_cal;
MissClassRate_cv = 100-Accuracy_cv;
Sen_cal = diag(cal_conf)./sum(cal_conf,2);
Spec_cal = (repmat(sum(cal_conf(:)),[size(cal_conf),1])-sum(cal_conf,2)-sum(cal_conf,1)+diag(cal_conf))./(repmat(sum(cal_conf(:)),[size(cal_conf),1]) - sum(cal_conf,2));
Sen_cv = diag(cv_conf)./sum(cv_conf,2);
Spec_cv = (repmat(sum(cv_conf(:)),[size(cv_conf),1])-sum(cv_conf,2)-sum(cv_conf,1)+diag(cv_conf))./(repmat(sum(cv_conf(:)),[size(cv_conf),1]) - sum(cv_conf,2));

Model.CM_stats.Accuracy_cal               = Accuracy_cal;
Model.CM_stats.Accuracy_cv                = Accuracy_cv;
Model.CM_stats.MissClassificationRate_cal = MissClassRate_cal;
Model.CM_stats.MissClassificationRate_cv  = MissClassRate_cv;
Model.CM_stats.Sensibility_cal            = Sen_cal;
Model.CM_stats.Sensibility_cv             = Sen_cv;
Model.CM_stats.Specificity_cal            = diag(Spec_cal);
Model.CM_stats.Specificity_cv             = diag(Spec_cv);

%% ROC curve
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave == 0
    Xroc = zeros(size(Model.Ycv(:,j),1)+1,size(u,1));
    Yroc = zeros(size(Model.Ycv(:,j),1)+1,size(u,1));
    AUC = zeros(size(u,1),1);
    for j=1:1:uu
        [Xroc(:,j),Yroc(:,j),~,AUC(j)] = perfcurve(Model.Yc(:,j),Model.Ycv(:,j),1);
    end
    
    Model.ROC.FPrate = Xroc;
    Model.ROC.TPrate = Yroc;
    Model.ROC.AUC    = AUC;
end
