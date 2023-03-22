function Xn=norma_DTLab(X,type)
% Xn=norma_DTLab(X,type);
% Does Normalization on X matrix
% Depending on the type value, a different normalization is proceeded
%
% USAGE :
% Xn=norma_DTLab(X,type);
%
% INPUT :
% X : Data matrix
% type : Type of normalization (0 = Length norm, 1 = Area norm, 2 = Max norm)
%
% OUTPUT :
% Xn : X matrix transformed by Normalization

switch type
    case 0 % Lenght
        n = vecnorm(X');
        Xn = X./n';
    case 1 % Lenght
        m = min(X,[],2);
        n = sum(X-m,2);
        Xn = X./n;
    case 2 % Lenght
        m = max(X,[],2);
        Xn = X./m;
end
