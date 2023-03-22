function dX = SG_Deriv_DTLab( X, winsize)
% dX = SG_Deriv_DTLab( X, winsize)
% Does Savitzky–Golay derivative (1st polynomial) on X matrix
%
% USAGE :
% dX = SG_Deriv_DTLab( X, winsize);
%
% INPUT :
% X : Data matrix
% winzise : Savitzky–Golay winsize filter
%
% OUTPUT :
% dX : X matrix transformed Savitzky–Golay derivative

 [Numer, Denom] = CoeffDer1_SG( winsize );
 
 [n, p] = size(X);
 halfstep=(winsize-1)/2;
 index = p - halfstep;
 
dX(:,1:halfstep) = zeros(n,halfstep);
for i = (halfstep + 1) : index
    S = Numer(halfstep+1) * X(:,i);
    for N = 1 : halfstep
        S = S + Numer(N) * (-X(:,i-N) + X(:,i+N));
    end
    dX(:,i) = (S / Denom);
 end
 dX(:,index+1:p) = zeros(n,halfstep);

function [ Numer, Denom ] = CoeffDer1_SG( winsize )

M=(winsize-1)/2;

Numer(M+1-(0:M)) = -1*(0:M);
Numer((M:2*M)+1) = (M:2*M)-M;

Denom = sum(Numer.*Numer);
