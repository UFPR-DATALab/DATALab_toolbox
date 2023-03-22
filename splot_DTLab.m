function [cov,corr,Ind] = splot_DTLab(X,t,options)
% Does S_Plot to visualise potentially interesting variables 
%
% USAGE :
% [corr,cov] = splot_DATALab(X,t,options)
%
% INPUT :
% X : Data matrix
% t : Corresponding scores eg.: model.loads{1,1}
% options.plot = display plot 'on' / 'off'
% options.treatment = 1 (mean centering) / 0 (none)
%
% OUTPUT :
% corr : CORR(t, X)
% cov : COV(t, X)
% Ind: Indexes of variables with n standard deviation values

%% Check options
if isnumeric(X) == 0
    X = X.data;
end

if exist('options','var')
    if isfield(options,'plot')
        plot = options.plot;
    else
        plot = 0;
    end
    if isfield(options,'treatment')
        ppx = options.treatment;
    else
        ppx = 1;
    end
else
    plot = 'on';
    ppx = 1;
end

%% Mean centering
if ppx == 1
    [nR,~]=size(X);
    mX = mean(X);
    Xpp = X - ones(nR,1)*mX;
else
    Xpp = X;
end

%% S-plot

N=size(Xpp,1);
LVs=size(t,2);
cov = zeros(LVs,size(Xpp,2));
corr = zeros(LVs,size(Xpp,2));

for i=1:LVs
    % Scale scores
    T = t(:,i)*inv(t(:,i)'*t(:,i));
    
    % Calculate the s-plot values.
    cov(i,:) = (T' * Xpp) ./ (N - 1);
    corr(i,:) = cov(i,:) ./ (std(T) .* std(Xpp));    
end

corr_vec = corr(1,:);
Point=[1:size(corr_vec,2)];
std_XrP=std(corr_vec);
std_XrN=-std(corr_vec);

% Calculate 1 to 3 x SDs
for j=1:3
    
    titiP=Point(corr_vec>j*std_XrP)';
    SD{j}=size(titiP,1);
    
    Ind{j}=titiP;
    Values{j}=corr_vec(1,titiP)';
    
    titiN=Point(corr_vec<j*std_XrN)';
    SD{j}=SD{j}+size(titiN,1);
    
    Ind{j}=[Ind{j};titiN];
    Values{j}=[Values{j};corr_vec(1,titiN)'];
    
end

if strcmp(plot,'on') == 1
    % Plot
    figure;
    hold on;
    title('S-plot: CORR(t, X) vs. COV(t, X)');
    xlabel('COV(t, X)');
    ylabel('CORR(t, X)');
    
    % Scatter plot the values
    htmp = scatter(cov(1,:), corr(1,:), 30, [0, 0, 0], 'filled');
       
    % Plot the s1 SD outliers in colour
    hold on;
    scatter(cov(1,Ind{1}),corr(1,Ind{1}),40,'g'); axis tight;
    scatter(cov(1,Ind{2}),corr(1,Ind{2}),40,'b'); axis tight;
    scatter(cov(1,Ind{3}),corr(1,Ind{3}),40,'r'); axis tight;
end

