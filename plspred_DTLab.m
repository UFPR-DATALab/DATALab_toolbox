function [Model_pred] = plspred_DTLab(Model,Xval,options,Yval)
% [Model_pred] = plspred_DTLab(Model,Xval,options,Yval);
%
% INPUT :
% Model : Output from plscal_DTLab
% Xval : Prediction/Validation X
% options.sampName : String matrix with sample names
% options.plots = display plot (0 (none) / 1 (some - DEFAULT) / 2 (all))
% Yval : Validation Y (For Blind test)
%
%
% OUTPUT:
% Model_pred
% Model_pred.Ypred : Predicted responses
% Model_pred.DModXn_pred: Normalized predicted DModX
% Model_pred.T2_pred: Predicted Hotteling
% Model_pred.RMSEV  = RMSEV for blind test
% Model_pred.R2V  = R²V for blind test
%
% version 1.0
% DATALab - Data Science in Chemistry and Chemometrics Laboratory
% 

%% Check options
if exist('options','var')
    if isfield(options,'plots')
        plots = options.plots;
    else
        plots = 1;
    end
else
    plots = 1;
end

%% Model data recovery
ppx = Model.ppx;
ppy = Model.ppy;
LV = Model.LV;
b = Model.beta;
Ppred = Model.Xloading;
XS = Model.Xscores;
Yc = Model.Yc;
W = Model.stats.W;
Tcrit = Model.stats.Tcrit;
Qcrit = Model.stats.Qcrit;
% s0 = Model.stats.s0;

%% Pre-treatment
if ppx == 1
    Xpp_v = meancentpred_DTLab(Xval,Model.mX);
elseif ppx == 2
    Xpp_v = autoscalpred_DTLab(Xval,Model.mX,Model.sX);
elseif ppx == 3
    Xpp_v = paretopred_DTLab(Xval,Model.mX,Model.sX);
else
    Xpp_v = Xval;
end


%% Prediction
Ypred = [ones(size(Xpp_v,1),1),Xpp_v]*b;

if ppy == 1
    Yp = Ypred + mean(Yc);
elseif ppy == 2
    Yp = Ypred.*std(Yc) + mean(Yc);
elseif ppy == 3
    Yp = Ypred.*sqrt(std(Yc)) + mean(Yc);
else
    Yp = Ypred;
end

if nargin > 3 %for blind test TO BE IMPLEMENTED
    rmsev(1) = sqrt(mean((Yval-Yp).^2));
    r = corrcoef(Yval,Yp);
    r2v = r.*r;
    r2v = r2v(1,2);
end
%% Diagnostic tools
Tpred = Xpp_v*W;
Xrec = Tpred*Ppred';

%----------------------------- Not necessary!!!------------------------------------------------------
% if ppx == 1
%     Xrec = Xrec + mean(Xval(:,Model.var_sel));
% elseif ppx == 2
%     Xrec = Xrec.*(ones(size(Xrec,1),1)*std(Xval(:,Model.var_sel))) + mean(Xval(:,Model.var_sel));
% elseif ppx == 3
%     Xrec = Xrec.*(ones(size(Xrec,1),1)*sqrt(std(Xval(:,Model.var_sel)))) + mean(Xval(:,Model.var_sel));
% end
%----------------------------------------------------------------------------------------------------

t2_pred = sum( bsxfun(@rdivide, abs(Tpred).^2, var(XS,[],1)) , 2);

ex = Xpp_v-Xrec;


% dx_pred = sqrt(sum(ex.^2,2)/(size(Xpp_v,2) - LV));
% dx_predn = dx_pred./s0;
for i=1:size(ex,1)
    Q_pred(i,1) = ex(i,:)*ex(i,:)';
end
qind = find(Q_pred>Qcrit);
tind = find(t2_pred>Tcrit);

% %% Report
Model_pred.Ypred = Yp;
Model_pred.Q_pred = Q_pred;
Model_pred.T2_pred  = t2_pred;
if nargin > 3
    Model_pred.RMSEV  = rmsev;
    Model_pred.R2V  = r2v;
end

if isfield(options,'sampName')
    Am = options.sampName;
else
    Am = 1:size(Xval,1);
    Am = cellstr(num2str(Am'));
end

%% Plots
if plots == 2
    %DModX figure
    axes1 = axes('Parent',figure);
    plot(Q_pred,'-o','MarkerFaceColor','auto'); hold on;
    if ~isempty(qind)
        scatter(qind,Q_pred(qind),40,'filled');hold on
                text(qind,Q_pred(qind),Am(qind), 'Interpreter', 'none');
    end
    line([0 size(Xval,1)],[Qcrit Qcrit],'Color','k','LineStyle','--');axis tight
    axesLimits = ylim(axes1)*1.2;
    ylim (axesLimits)
    xlabel 'Predicted samples'
    ylabel 'Normalized distance to the model in X'
    title 'Residuals prediction'

    %Hoteling figure
    axes1 = axes('Parent',figure);
    plot(t2_pred,'-or','MarkerFaceColor','auto'); hold on;
    if ~isempty(tind)
        scatter(tind,t2_pred(tind),40,'b','filled');hold on
        text(tind,t2_pred(tind),Am(tind), 'Interpreter', 'none');
    end
    line([0 size(Xval,1)],[Tcrit Tcrit],'Color','k','LineStyle','--');axis tight
    axesLimits = ylim(axes1)*1.2;
    ylim (axesLimits)
    xlabel 'Predicted samples'
    ylabel 'Predicted T^2'
    title 'Hotteling T^2 prediction'
end

if plots >= 1
  if nargin > 3
    for y = 1:1:size(Yval,2)
        %Validation figure
        axes1 = axes('Parent',figure);
        sc1=scatter(Model.Yc(:,y),Model.Yfit(:,y),40,'filled'); hold on;
        sc2=scatter(Yval(:,y),Yp(:,y),40,'r','filled');
        axesLimits1 = xlim(axes1);
        ax = axis;
        xplot1 = linspace(ax(1),ax(2));
        fitResults1 = polyfit(Yval(:,y),Yp(:,y), 1);
        yplot1 = polyval(fitResults1, xplot1); hold on;
        fitLine1 = plot(xplot1,yplot1,'DisplayName','   linear','Parent',axes1,...
            'Tag','linear',...
            'Color',[0.47 0.67 0.19],'linestyle','--');
        ax = axis;
        [~,im] = max(ax);
        Rfl = line([ax(im-1) ax(im)],[ax(im-1) ax(im)],'Color',[0.91 0.41 0.10],'linestyle','-');
%         Rfl = refline(1,0);
%         Rfl.Color = [0.91 0.41 0.10];
%         Rfl.LineStyle = '-';
        
        xlabel 'Reference Y'
        ylabel 'Predicted Y'
        title 'Validation correlation'
        annotation('textbox',...
            [0.16 0.73 0.3 0.15],...
            'String',{['RMSEV = ' num2str(rmsev(1,y))],['R² = ' num2str(r2v(1,y))]},...
            'FitBoxToText','off',...
            'EdgeColor','none');
%         legend('Calibration Samples', 'Validation Samples','Adjusted Line','Reference Line')
        box 'on'
        grid 'on'
    end
    end
    
    %DModX vs T2 figure
    grp = ones(size(Xval,1),1);
    int = intersect(qind,tind);    
    if ~isempty(qind)
        grp(qind) = 2;
        if ~isempty(tind)
            grp(tind) = 3;
            if ~isempty(int)
                grp(int) = 4;
            end
        end
    end
    u = unique (grp);
    if ~isempty(u(u==4)) && ~isempty(u(u==3)) && ~isempty(u(u==2))
        leg = {'Predicted samples', 'Higher DModx', 'Higher T^2', 'Higher DModX and T^2'};
    elseif ~isempty(u(u==4)) && ~isempty(u(u==3)) && isempty(u(u==2))
        leg = {'Predicted samples', 'Higher T^2', 'Higher DModX and T^2'};
    elseif ~isempty(u(u==4)) && isempty(u(u==3)) && ~isempty(u(u==2))
        leg = {'Predicted samples', 'Higher DModx', 'Higher DModX and T^2'};
    elseif ~isempty(u(u==4)) && isempty(u(u==3)) && isempty(u(u==2))
        leg = {'Predicted samples', 'Higher DModX and T^2'};
    elseif isempty(u(u==4)) && ~isempty(u(u==3)) && ~isempty(u(u==2))
        leg = {'Predicted samples', 'Higher DModx', 'Higher T^2'};
    elseif isempty(u(u==4)) && isempty(u(u==3)) && ~isempty(u(u==2))
        leg = {'Predicted samples', 'Higher DModx'};
    elseif isempty(u(u==4)) && ~isempty(u(u==3)) && isempty(u(u==2))
        leg = {'Predicted samples', 'Higher T^2'};
    else
        leg = {'Predicted samples'};
    end
    
    figure;
    scatter(Q_pred(grp==1),t2_pred(grp==1),25,'r','filled');hold on
    scatter(Q_pred(grp==2),t2_pred(grp==2),25,'g','filled');
    scatter(Q_pred(grp==3),t2_pred(grp==3),25,'b','filled');
    scatter(Q_pred(grp==4),t2_pred(grp==4),25,[.5 0 .5],'filled');
%     gscatter(Q_pred,t2_pred,grp);
    minx = min(Q_pred);
    maxx = max(Q_pred);
    miny = min(t2_pred);
    maxy = max(t2_pred);
    xlim([min([minx*1.2 0]) maxx*1.2]);
    ylim([min([miny*1.2 0]) maxy*1.2]);
    if ~isempty(qind)
        text(Q_pred(qind),t2_pred(qind),Am(qind), 'Interpreter', 'none');
        if ~isempty(tind)
            text(Q_pred(tind),t2_pred(tind),Am(tind), 'Interpreter', 'none');
            if ~isempty(int)
                text(Q_pred(int),t2_pred(int),Am(int), 'Interpreter', 'none');
            end
        end
    end
    if Qcrit<max(maxx)
        line([Qcrit Qcrit],[min([miny*1.2 0]) maxy*1.2],'Color','k','LineStyle','--');
        leg = [leg, {'D_c_r_i_t'}];
    end
    if Tcrit<max(maxy)
        ax = axis;
        tline = line([ax(1) ax(2)],[Tcrit Tcrit],'Color','k','linestyle','--');
%         tline = refline(0,Tcrit*2.5);
%         tline.Color = 'k';
%         tline.LineStyle = '--';
        leg = [leg, {'T²_c_r_i_t'}];
    end
    xlabel 'DModXn prediction'
    ylabel 'Predicted T^2'
    title 'DModXn v. Hotteling T^2 prediction'
    box on
    grid on
    legend (leg)
end
