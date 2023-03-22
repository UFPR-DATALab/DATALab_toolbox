function [Model] = PCA_DTLab(X,options)
% [Model] = PCA_DTLab(X,options);
%
% INPUT :
% X : Data matrix
% options.sampName : cell containing Sample Names
% options.varName : cell containing Variable Names
% options.grpClass : cell containing Group Classes
%
% OUTPUT:
% Model.P = Loadings;
% Model.T = Scores;
% Model.VarExp = Explained Variance;
% Model.t2 = Hotelling's T²;
% Model.Q = Residuals;
%
% USAGE:
% pca - Statistics and Machine Learning Toolbox
%
% version 1.0
% DATALab - Data Science in Chemistry and Chemometrics Laboratory
% 

%% PCA
[n, p] = size(X);
% Number of components
[~,s,~] = svd(X,'econ');
explained = (diag(s.^2)/sum(diag(s.^2)))*100;
figure;
plot(explained,'ob:')
xlabel 'Number of Components'
ylabel 'Explained Variance (%)'

lv_txt = ('How many components? ');
ncomp = input(lv_txt);

% Calculate PCA with the corrected number of PCs
% [coeff,score,latent,tsquared,explained] = pca(X);
[u,s,v] = svd(X,'econ');
Coeff = v;
Scores = u*s;
explained = (diag(s.^2)/sum(diag(s.^2)))*100;
q = X - Scores(:,1:ncomp)*Coeff(:,1:ncomp)';
for i=1:n
    Q(i) = q(i,:)*q(i,:)';
end

latent = diag(s.^2)./n;
% latent = latent(1:ncomp);
q = sum(latent > max(n,p)*eps(latent(1)));
standScores = bsxfun(@times, Scores(:,1:q), 1./sqrt(latent(1:q,:))');
tsquared = sum(standScores.^2, 2);

% [Coeff,Scores,~,tsquared,explained]=pca(X,'NumComponents',ncomp,'Centered',false);

% Critical Values
Fcrit = finv(0.95,ncomp,size(X,1)-ncomp);
Tcrit = (ncomp*(size(X,1)-1)/(size(X,1)-ncomp))*Fcrit;
dx = sqrt(sum(Q.^2,2)/(size(X,2) - ncomp));

% Plot Scores
figure;
if exist('options','var')
    if isfield(options,'grpClass')
        [C,~,IC] = unique(options.grpClass);
        for i=1:size(C,1)
            hold on;
            scatter(Scores(IC==i,1),Scores(IC==i,2),25,'filled');
        end
    else
        scatter(Scores(:,1),Scores(:,2),25,'filled');
    end
else
    scatter(Scores(:,1),Scores(:,2),25,'filled');
end
%xlabel(['PC1 ' '(' num2str(round(explained(1,1),2)) '%' ')']); 
%ylabel(['PC2 ' '(' num2str(round(explained(2,1),2)) '%' ')']);
xlabel(['PC1 ' '(' num2str(explained(1,1),2) '%' ')']); 
ylabel(['PC2 ' '(' num2str(explained(2,1),2) '%' ')']);
grid on;
box on;
ax = axis;
if ax(3)<0 && ax(4)>0
    line([ax(1),ax(2)], [0 0],'LineStyle','-');
end
if ax(1)<0 && ax(2)>0
    line([0 0],[ax(3),ax(4)],'LineStyle','-');
end
if exist('options','var')
    if isfield(options,'sampName')
        samp_labels = options.sampName;
        [m,n] = size(samp_labels);
        if ~iscellstr(samp_labels)
            samp_labels = cellstr(samp_labels);
        end
        if m>n
            for k=1:1:length(samp_labels)
                text(Scores(k,1),Scores(k,2)-0.04,samp_labels{k,1});
            end
        else
            for k=1:1:length(samp_labels)
                text(Scores(k,1),Scores(k,2)-0.04,samp_labels{1,k});
            end
        end
    end
end
if exist('options','var')
    if isfield(options,'grpClass')
        legend(C)
    end
end
axis tight

% Plot Loadings
figure;
scatter(Coeff(:,1),Coeff(:,2),25,[1 0 0],'filled');
xlabel(['PC1 ' '(' num2str(explained(1,1),2) '%' ')']); 
ylabel(['PC2 ' '(' num2str(explained(2,1),2) '%' ')']);
grid on;
box on;
ax = axis;
if ax(3)<0 && ax(4)>0
    line([ax(1),ax(2)], [0 0],'LineStyle','-');
end
if ax(1)<0 && ax(2)>0
    line([0 0],[ax(3),ax(4)],'LineStyle','-');
end
if exist('options','var')
    if isfield(options,'varName')
        var_labels = options.varName;
        [m,n] = size(var_labels);
        if ~iscellstr(var_labels)
            if ischar(var_labels)
                var_labels = cellstr(var_labels);
            elseif isnumeric(var_labels)
                if m>n
                    var_labels = num2str(var_labels);
                else
                    var_labels = num2str(var_labels');
                end
                var_labels = cellstr(var_labels);
            end
        end
        [m,n] = size(var_labels);
        if m>n
            for k=1:1:length(var_labels)
                text(Coeff(k,1),Coeff(k,2)-0.04,var_labels{k,1});
            end
        else
            for k=1:1:length(var_labels)
                text(Coeff(k,1),Coeff(k,2)-0.04,var_labels{1,k});
            end
        end
    end
end
axis tight

% Plot Bi-plot
figure;

% Test to see if scores and loadings are adjusted------------------------
% Force each column of the coefficients to have a positive largest element.
% This tends to put the large var vectors in the top and right halves of
% the plot.
[p,d] = size(Coeff);
[~,maxind] = max(abs(Coeff),[],1);
colsign = sign(Coeff(maxind + (0:p:(d-1)*p)));
Coeff_bplt = bsxfun(@times,Coeff,colsign);


maxCoefLen = sqrt(max(sum(Coeff_bplt(:,1:2).^2,2)));
Scores_bplt = bsxfun(@times, maxCoefLen.*(Scores ./ max(abs(Scores(:)))), colsign);
%----------------------------------------------------------------    

% scatter(Scores(:,1)./norm(Scores(:,1)),Scores(:,2)./norm(Scores(:,2)),25,'filled');
scatter(Scores_bplt(:,1),Scores_bplt(:,2),25,'filled');
hold on
% scatter(Coeff(:,1),Coeff(:,2),25,[1 0 0],'filled')
scatter(Coeff_bplt(:,1),Coeff_bplt(:,2),25,[1 0 0],'filled')
xlabel(['PC1 ' '(' num2str(explained(1,1),2) '%' ')']); 
ylabel(['PC2 ' '(' num2str(explained(2,1),2) '%' ')']);
grid on;
box on;
ax = axis;
if ax(3)<0 && ax(4)>0
    line([ax(1),ax(2)], [0 0],'LineStyle','-');
end
if ax(1)<0 && ax(2)>0
    line([0 0],[ax(3),ax(4)],'LineStyle','-');
end
if exist('options','var')
    if isfield(options,'sampName')
        samp_labels = options.sampName;
        [m,n] = size(samp_labels);
        if ~iscellstr(samp_labels)
          samp_labels = cellstr(samp_labels);
        end 
        if m>n        
            for k=1:1:length(samp_labels)
              text(Scores(k,1),Scores(k,2)-0.04,samp_labels{k,1});
            end
          else
            for k=1:1:length(samp_labels)
              text(Scores(k,1),Scores(k,2)-0.04,samp_labels{1,k});
            end
        end
    end
    if isfield(options,'varName')
        var_labels = options.varName;
        [m,n] = size(var_labels);
        if ~iscellstr(var_labels)
            if ischar(var_labels)
                var_labels = cellstr(var_labels);
            elseif isnumeric(var_labels)
                if m>n
                    var_labels = num2str(var_labels);
                else
                    var_labels = num2str(var_labels');
                end
                var_labels = cellstr(var_labels);
            end
        end
        [m,n] = size(var_labels);
        if m>n
            for k=1:1:length(var_labels)
                text(Coeff(k,1),Coeff(k,2)-0.04,var_labels{k,1});
            end
        else
            for k=1:1:length(var_labels)
                text(Coeff(k,1),Coeff(k,2)-0.04,var_labels{1,k});
            end
        end
    end
end
axis tight

% Plot Residuals
figure;
Res_grp = cell(size(X,1),1);
Res_grp(tsquared>Tcrit,1) = {'High T²'};
Res_grp(Q'>dx,1) = {'High Residuals'};
Res_grp(tsquared>Tcrit & Q'>dx,1) = {'High T² and High Residuals'};
Res_grp(tsquared<Tcrit & Q'<dx,1) = {'Samples'};
[C,~,IC] = unique(Res_grp);
for i=1:size(C,1)
    hold on;
    scatter(tsquared(IC==i),Q(IC==i),25,'filled');
end
xlabel('Hotelling T²'); 
ylabel('Residuals');
grid on;
box on;
if exist('options','var')
    if isfield(options,'sampName')
        samp_labels = options.sampName;
        [m,n] = size(samp_labels);
        if ~iscellstr(samp_labels)
          samp_labels = cellstr(samp_labels);
        end 
        if m>n        
            for k=1:1:length(samp_labels)
%                if ~contains(Res_grp(k),'Samples')
                if ~ismember('Samples',Res_grp(k))
                    text(tsquared(k),Q(k)-0.04,samp_labels{k,1});
                end
            end
          else
            for k=1:1:length(samp_labels)
%          if ~contains(Res_grp(k),'Samples')
                if ~ismember('Samples',Res_grp(k))
                    text(tsquared(k),Q(k)-0.04,samp_labels{1,k});
                end
            end
        end
    end
end
ax = axis;
line([ax(1),ax(2)], [dx dx],'LineStyle','--','Color','w');
line([Tcrit Tcrit],[ax(3),ax(4)],'LineStyle','--','Color','w');
ax = axis;
line([ax(1),ax(2)], [dx dx],'LineStyle','--','Color','k');
line([Tcrit Tcrit],[ax(3),ax(4)],'LineStyle','--','Color','k');
legend(C)
axis tight

% Save data
Model.P = Coeff;
Model.T = Scores;
Model.VarExp = explained;
Model.t2 = tsquared;
Model.Q = Q;

