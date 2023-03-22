function [trshld, prob] = plsda_bayseatnrshld_DTLab(Y,Ypred,options)
% [trshld, prob] = plsda_bayseatnrshld_DTLab(Y,Ypred,options)
%
% Baysean approach to determine the PLS-DA class limits



class = unique(Y);

llim = min(Ypred);
hlim = max(Ypred);
prec   = abs(hlim-llim)/100; % Precision of gaussian curve
if prec == 0; prec = 1; end
inrng  = (llim:prec:hlim)'; % Range vector of gaussian curve

% Calculate distribution probabilities
for i=1:length(class)
    s =  std(Ypred(Y==class(i)));
    m = mean(Ypred(Y==class(i)));
    p(:,i) = 1./(sqrt(2*pi)*(s./prec)) * exp(-0.5*(([0:(length(inrng)-1)]'-((m-llim)./prec))./(s./prec)).^2);
end

% Normalization step
p = p + 1.e-11;  % to prevent norm = 0.
pn = sum(p,2);
pn(pn==0) = inf;
prob = p ./ (pn*ones(1,size(p,2)));

prob(~isfinite(prob)) = 0;  %default probability if not finite
[~, pos] = min(abs(diff(prob')),[],2);
% [~,pos] = findpeaks(abs(diff(prob'))*-1);
if size(pos,2)>1
    pos_corr = pos(1);
    if inrng(pos_corr)<0
        pos_corr = pos(2);
    end
else
    pos_corr = pos;
end

prob = [inrng(:) prob];     %tag range of y values on as first column

% Calculate Threshold
trshld = inrng(pos_corr);

if options.plots == 1
    figure;
    for i=1:length(class)
        [x(:,i),n(:,i)] = hist(Ypred(Y==class(i)));
    end
    bar(n,x,200);hold on;

    xlabel 'PLS-Predicted Value'
    ylabel 'Number of Samples'
    ax = axis;
    line([trshld trshld],[ax(3) ax(4)],'Color','k','LineStyle','--');
    plot(prob(:,1),prob(:,2:end)*ax(4))
end
    


