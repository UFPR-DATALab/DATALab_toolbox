function Xali=SG_ali_DTLab(X,winsize)
% [Xali]=SG_ali_DTLab(X,winsize)
% Does Savitzky–Golay smoothing on X matrix
%
% USAGE :
% [Xali]=SG_ali_DTLab(X,winsize)
%
% INPUT :
% X : Data matrix
% winzise : Savitzky–Golay winsize filter
%
% OUTPUT :
% Xali : smoothed X matrix 

[ Numer, Denom ] = CoeffSmooth_SG( winsize );
 
[n, p] = size(X);
 
demipas=(winsize-1)/2;
indice = p - demipas;
 
    % First points unchanged
	Xali(:,1:demipas) = X(:,1:demipas);
 
    for i = (demipas + 1) : indice
        Somme = Numer(demipas+1) * X(:,i);
        for N = 1 : demipas
            Somme = Somme + Numer(N) * (X(:,i-N) + X(:,i+N));
        end
        Xali(:,i) = (Somme / Denom);
 
    end
 
    % Last points unchanged
    Xali(:,indice:p) = X(:,indice:p);

function [ Numer, Denom ] = CoeffSmooth_SG( winsize )
 
M=(winsize -1)/2;
 
for i = 0:M
        Numer(M-i+1) =  (3 * M * M + 3 * M - 1 - 5 * (i) * (i));
end
    
for i = M:2*M
        Numer(i+1) =  (3 * M * M + 3 * M - 1 - 5 * (i-M) * (i-M));
end
 
Denom = (2 * M + 3) * (2 * M + 1) * (2 * M - 1)/3;
