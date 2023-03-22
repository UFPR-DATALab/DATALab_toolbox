function CM = confusionmat(g,ghat,varargin)
%% Help - Need to do later

C = unique(g,'rows');
[ac,~] = size(C);
Chat = unique(ghat,'rows');
[achat,~] = size(Chat);

CM = zeros(ac,achat);

for i=1:ac
    A = g == C(i,1);
    for j=1:ac
        CM(i,j) = sum(ghat(A) == C(j,1));
    end
    if sum(ghat(A) == 0)~=0
        CM(i,achat) = CM(i,achat)+sum(ghat(A) == 0);
    end
end
CM = CM';
