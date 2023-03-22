function Xpar=pareto_DTLab(X)
% Xpar=pareto_DTLab(X);
% Does Pareto scaling on X matrix
% Should I make without mean center(???)
%
% USAGE :
% Xpar=par_DTLab(X);
%
% INPUT :
% X : Data matrix
%
% OUTPUT :
% Xpar : X matrix transformed by Pareto scaling

[n,m]=size(X);

% Mean center
mX = mean(X);
Xmc = X - ones(n,1) * mX;

% Pareto
% Xpar= Xmc ./ (ones(n,1) * (sqrt(std(Xmc))));
Xpar = zeros(n,m);
z = (std(Xmc)==0);
d = (ones(n,1) * (sqrt(std(Xmc))));
Xpar(:,~z)= Xmc(:,~z) ./ d(:,~z);
