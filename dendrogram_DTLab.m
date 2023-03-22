function [xtl] = dendrogram_DTLab(X, Z, cutoff)
% [xtl] = dendrogram_DTLab(X, Z, cutoff);
%
% INPUT :
% X : Data matrix
% Z : Linkage matrix obtained from HCA_DTLab
% cutoff : cut position in dendograms (default 0.5*max(linkage matrix))
%
% OUTPUT:
% xtl = X tick label from dendrogram
% version 1.0
% DATALab - Data Science in Chemistry and Chemometrics Laboratory
% 

%% Calculate tree
[Roots, Clusters, Number] = Analysis_tree(Z,cutoff);

%% Create figure
% Axes
h = figure;
ax1 = axes;
set(ax1,'Ygrid','on','TickLength',[0 0],'YMinorGrid','on');
title('Hierarchical clustering & dendrogram',...
      'Fontsize',   12);
ylabel(ax1,'Dissimilarity');
hold(ax1,'on');

% Number of elements
N = max(max(Z(:,1:2)))/2+1;

% X width
M = size(X,2);

% Number of roots
R = numel(Roots);

% Number of Z or nodes
T = size(Z,1);

% Roots coordinates
AbscissaRoots  = 1:R;
OrdinatesRoots = zeros(1,R);

% Nodes coordinates
AbscissaNodes  = zeros(1,N);
OrdinatesNodes = zeros(1,N);

% Axes limits
Xlim = [0 N+1];
k = 5;
Ylim = [0 ceil(max(Z(:,3))*k)/k];
set(ax1,...
    'Xlim',     Xlim,...
    'Ylim',     Ylim,...
    'Fontsize', 9);

% Indices
Indices = 1:N;

% Colormap
Colors = colormap('jet');

% Extrema
Min = min(arrayfun(@(n)norm(X(n,:)),1:size(X,1)));
Max = max(arrayfun(@(n)norm(X(n,:)),1:size(X,1)));

% Y axes limit
Ylim = get(ax1,'Ylim');

% Clustering limit
Limit = mean([Z(end-Number+2,3) Z(end-Number+1,3)]);
plot([0.5 N+0.5],Limit*[1 1],'r--');

% Figure
Figure = get(ax1,'Parent');

%% Loop to create dendrogram
for t = 1:T
    
    % Nodes or members of the current cluster
    Node1 = Z(t,1);
    Node2 = Z(t,2);
    
    % First member of the current cluster
    if Node1 <= N
        
        % Index of the current node
        i = find(Roots==Node1);
        
        % Indices permutation
        Index = Indices(i);
        Indices(Indices==Node1) = Index;
        Indices(i) = Node1;
        
        % Coordinates of the current root
        x1 = AbscissaRoots(i);
        y1 = OrdinatesRoots(i);
        
    else
        
        % Coordinates of the current node
        x1 = AbscissaNodes(Node1-N);
        y1 = OrdinatesNodes(Node1-N);
        
    end
    
    % Second member of the cluster
    if Node2 <= N
        
        % Index of the current node
        i = find(Roots==Node2);
        
        % Indices permutation
        Index = Indices(i);
        Indices(Indices==Node2) = Index;
        Indices(i) = Node2;
        
        % Coordinates of the current root
        x2 = AbscissaRoots(i);
        y2 = OrdinatesRoots(i);
        
    else
        
        % Coordinates of the current node
        x2 = AbscissaNodes(Node2-N);
        y2 = OrdinatesNodes(Node2-N);
        
    end
       
    % Coordinates of the current node
    y = Z(t,3);
    AbscissaNodes(t) = mean([x1 x2]);
    OrdinatesNodes(t) = y;
    
    % Data label
    Xticklabel = arrayfun(@(n)sprintf('%u',Indices(n)),1:N,'UniformOutput',false);
    set(ax1,'Xtick',1:N,'XtickLabel',Xticklabel);
    
    
    % Creation of a new graphical node
    plot(ax1,[x1 x1 x2 x2],[y1 y y y2],'b');
    set(ax1,'Xlim',Xlim,'Ylim',Ylim);
    
    drawnow();
end
xtl = Xticklabel;


