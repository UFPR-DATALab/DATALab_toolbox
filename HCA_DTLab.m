function Z = HCA_DTLab(X,option)
% Z = HCA_DTLab(X,options);
%
% INPUT :
% X : Data matrix
% options.dist : distance method ('euclidean', 'cityblock', 'mahalanobis' default = 'euclidean')
% options.link : linkage method ('average', 'single', 'complete', 'centroid', 'ward' default = 'ward')
% options.cutoff : cut position in dendograms (default 0.5*max(linkage matrix))
% options.sampName : cell containing Sample Names
%
% OUTPUT:
% Z : Linkage matrix
%
% USAGE:
% linkage - Statistics and Machine Learning Toolbox
% dendrogram_DTLab - DATALab toolbox
% Analysis_tree - DATALab toolbox
%
% version 1.0
% DATALab - Data Science in Chemistry and Chemometrics Laboratory
% 



%% HCA
% Options
if exist('options','var')
    if isfield(options,'dist')
        dist = options.dist;
    else
        dist = 'euclidean';
    end
    if isfield(options,'link')
        link = options.link;
    else
        link = 'ward';
    end
end
    
% Linkage step
Z = linkage(X,link,dist);
if exist('options','var')
    if isfield(options,'cutoff')
        cutoff = options.cutoff;
    else
        cutoff = 0.5*max(Z(:,3));
    end
end

% Dendogram
[xtl] = dendrogram_DTLab(X, Z, cutoff);
if exist('options','var')
    if isfield(options,'sampName')
        sampName = options.sampName;
    else
        supName = 1:size(X,1);
        sampName = cellstr(num2str(supName'));
    end
end
for i=1:size(xtl,2)
    xtl{1,i} = str2num(xtl{1,i});
end
xtl = cell2mat(xtl);
if iscell(sampName)
    sampName = cellstr(sampName);
end
xticklabels(sampName(xtl))
% xtickangle(90)
