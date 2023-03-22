function [Model] = plscal_DTLab(Xcal,Ycal,LV,options)
% [Model] = plscal_DTLab(Xcal,Ycal,LV,options);
%
% INPUT :
% Xcal : Calibration X
% Ycal : Calibration Y
% LV : Maximum latent variables to perform the model
% options.sampName = Cellstring matrix with sample names
% options.plots = display plot (0 (none) / 1 (some - DEFAULT) / 2 (all))
% options.ppx = Pre-treatment on X (0 (none - default) / 1 (Mean Center) / 2 (Autoscaling) / 3 (Pareto))
% options.ppy = Pre-treatment on Y (0 (none - default) / 1 (Mean Center) / 2 (Autoscaling) / 3 (Pareto))
% options.cv = {'vet',k};
% options.pls_type = type of PLS for more than one reponse (1 (PLS1) / 2
% (PLS2))
%
% OUTPUT:
% Model
% Model.beta : Regression vector
% Model.Xloading : Loadings of X resulted from PLS
% Model.LV : Selected latent variables for prediction
% Model.ppx : Pre-treatment of Xcal
% Model.ppy : Pre-treatment of Ycal
% Model.stats : Statistics parameters of the model (DModX, T²...)
%
% USAGE:
% plsregress - Statistics and Machine Learning Toolbox
%
% version 1.0
% DATALab - Data Science in Chemistry and Chemometrics Laboratory
% 

%% Check options
h1 = waitbar(0,'PLS calculation:');
waitbar(1/9,h1,['PLS calculation: ', sprintf('%0.1u',1) '/', sprintf('%0.1u',9)]);

if exist('options','var')
    if isfield(options,'cv')
        cv = options.cv{1,1};
        if strcmp(options.cv{1,1},'vet') == 1
            cv = 'cv';
            kfold = options.cv{1,2};
        end
    else
        cv = 'cv';
        kfold = 7;
    end
    if isfield(options,'plots')
        plots = options.plots;
    else
        plots = 0;
    end
    if isfield(options,'ppx')
        ppx = options.ppx;
    else
        ppx = 1;
    end
    if isfield(options,'ppy')
        ppy = options.ppy;
    else
        ppy = 1;
    end
    if isfield(options,'sampName')
        Am = options.sampName;
    else
        Am = 1:size(Xcal,1);
        Am = cellstr(num2str(Am'));
    end
    if isfield(options,'pls_type')
        ptype = options.pls_type;
    else
        ptype = 1;
    end
else
    cv = 'cv';
    kfold = 7;
    plots = 1;
    ppx = 0;
    ppy = 0;
    Am = 1:size(Xcal,1);
    Am = cellstr(num2str(Am'));
    ptype = 1;
end

% Am = options.sampName;

%% Pre-treatment
waitbar(2/9,h1,['PLS calculation: ', sprintf('%0.1u',2) '/', sprintf('%0.1u',9)]);

if ppx == 1
    Xpp = meancent_DTLab(Xcal);
elseif ppx == 2
    Xpp = autoscal_DTLab(Xcal);
elseif ppx == 3
    Xpp = pareto_DTLab(Xcal);
else
    Xpp = Xcal;
end

if ppy == 1
    Ypp = meancent_DTLab(Ycal);
elseif ppy == 2
    Ypp = autoscal_DTLab(Ycal);
elseif ppy == 3
    Ypp = pareto_DTLab(Ycal);
else
    Ypp = Ycal;
end

%% First calibration
waitbar(3/9,h1,['PLS calculation: ', sprintf('%0.1u',3) '/', sprintf('%0.1u',9)]);

% Calibration
try
    if ptype == 1
    [~,~,~,~,~,~,MSE] = plsregress(Xpp,Ypp,LV);
    else
        for i=1:size(Ypp,2)
            [~,~,~,~,~,~,MSE{i}] = plsregress(Xpp,Ypp(:,i),LV);
        end
    end
catch
    LV = 10;
    if ptype == 1
    [~,~,~,~,~,~,MSE] = plsregress(Xpp,Ypp,LV);
    else
        for i=1:size(Ypp,2)
            [~,~,~,~,~,~,MSE(:,i)] = plsregress(Xpp,Ypp(:,i),LV);
        end
    end
end
if ptype == 1
    rmse_lv(1,:) = MSE(2,2:end).^(1/2);
    [~,~,~,~,~,~,MSE] = plsregress(Xpp,Ypp,LV,cv,kfold);
    rmse_lv(2,:) = MSE(2,2:end).^(1/2);
        
    % LV selection
    figure;
    plot(rmse_lv',':o'); axis tight
    xlabel 'LVs'
    ylabel 'RMSE'
    legend 'RMSEC' 'RMSECV'
    
    lv_txt = ('How many latent variables? ');
    ncomp = input(lv_txt);
else
    for i=1:size(Ypp,2)
        rmse_lv(1+2*(i-1),:) = MSE{i}(2,2:end).^(1/2);
        [~,~,~,~,~,~,MSE{i}] = plsregress(Xpp,Ypp(:,i),LV,cv,kfold);
        rmse_lv(2+2*(i-1),:) = MSE{i}(2,2:end).^(1/2);
    end
    
    % LV selection
    for i=1:size(Ypp,2)
        figure;
        plot(rmse_lv([1+2*(i-1) 2+2*(i-1)],:)',':o'); axis tight
        xlabel 'LVs'
        ylabel 'RMSE'
        legend 'RMSEC' 'RMSECV'
        
        lv_txt = ('How many latent variables? ');
        ncomp(i) = input(lv_txt);
        close all
    end
end

%% Final Model
waitbar(4/9,h1,['PLS calculation: ', sprintf('%0.1u',5) '/', sprintf('%0.1u',9)]);

%--------------------------------------------------------------------------
% Does it needs this???
% [~,~,~,~,~,~,MSE,stats] = plsregress(Xpp,Ypp,LV);
% rmse(1,:) = MSE(2,2:end).^1/2;
% [~,~,~,~,~,~,MSE,stats] = plsregress(Xpp,Ypp,LV,cv,kfold);
% rmse(2,:) = MSE(2,2:end).^1/2;
% 
% % LV selection
% figure;
% plot(rmse',':o'); axis tight
% xlabel 'LVs'
% ylabel 'RMSE'
% legend 'RMSEC' 'RMSECV'
% 
% lv_txt = ('How many latent variables ?');
% ncomp = input(lv_txt);
%--------------------------------------------------------------------------

if ptype == 1 % PLS1
    [XL,~,XS,~,BETA,~,~,stats_final] = plsregress(Xpp,Ypp,ncomp);
    if ppy == 1
        Yc = (Ypp - stats_final.Yresiduals) + mean(Ycal(:));
    elseif ppy == 2
        Yc = (Ypp - stats_final.Yresiduals)*std(Ycal(:)) + mean(Ycal(:));
    elseif ppy == 3
        Yc = (Ypp - stats_final.Yresiduals)*sqrt(std(Ycal(:))) + mean(Ycal(:));
    else
        Yc = (Ypp - stats_final.Yresiduals);
    end
    rmse = zeros(2,size(Ycal,2));
    R = zeros(2,size(Ycal,2));
    for y=1:1:size(Ycal,2)
        rmse(1,y) = sqrt(mean((Ycal(:,y)-Yc(:,y)).^2));
        r = (corrcoef(Ycal(:,y),Yc(:,y))).^2;
        R(1,y) = r(1,2);
    end
else % PLS2
    rmse = zeros(2,size(Ycal,2));
    R = zeros(2,size(Ycal,2));
    for i = 1:size(Ypp,2)
        [XL{i},~,XS{i},~,BETA(:,i),~,~,stats_final{i}] = plsregress(Xpp,Ypp(:,i),ncomp(i));
        if ppy == 1
            Yc(:,i) = (Ypp(:,i) - stats_final{i}.Yresiduals) + mean(Ycal(:,i));
        elseif ppy == 2
            Yc(:,i) = (Ypp(:,i) - stats_final{i}.Yresiduals)*std(Ycal(:,i)) + mean(Ycal(:,i));
        elseif ppy == 3
            Yc(:,i) = (Ypp(:,i) - stats_final{i}.Yresiduals)*sqrt(std(Ycal(:,i))) + mean(Ycal(:,i));
        else
            Yc(:,i) = (Ypp(:,i) - stats_final{i}.Yresiduals);
        end        
        rmse(1,i) = sqrt(mean((Ycal(:,i)-Yc(:,i)).^2));
        r = (corrcoef(Ycal(:,i),Yc(:,i))).^2;
        R(1,i) = r(1,2);
    end
end

%% Perform CV - Venetian blind
waitbar(5/9,h1,['PLS calculation: ', sprintf('%0.1u',6) '/', sprintf('%0.1u',9)]);

if ptype == 1 % PLS1
    
    Ycv = zeros(size(Ypp,1),size(Ypp,2));
    
    for i=1:kfold
        %     (round(size(Ypp,1)/kfold))
        tst = (1+(i-1):kfold:size(Ypp,1));
        Xtstc = Xpp;
        Xtstc(tst,:) = [];
        Xtstp = Xpp(tst,:);
        Ytstc = Ypp;
        Ytstc(tst,:) = [];
        [~,~,~,~,b,~,~,~] = plsregress(Xtstc,Ytstc,ncomp);
        Ytstp = [ones(size(Xtstp,1),1),Xtstp]*b;
        Ycv(tst,:) = Ytstp;
    end
    if ppy == 1
        Ycv = Ycv + mean(Ycal(:,:));
    elseif ppy == 2
        Ycv = Ycv.*std(Ycal(:,:)) + mean(Ycal(:,:));
    elseif ppy == 3
        Ycv = Ycv.*sqrt(std(Ycal(:,:))) + mean(Ycal(:,:));
    end
    
    Rcv = zeros(1,y);
    for y=1:1:size(Ycal,2)
        rmse(2,y) = sqrt(mean((Ycal(:,y)-Ycv(:,y)).^2));
        r = (corrcoef(Ycal(:,y),Ycv(:,y))).^2;
        Rcv(1,y) = r(1,2);
    end
    R(2,:) = Rcv;

else % PLS2

    Ycv = zeros(size(Ypp,1),size(Ypp,2));
    Rcv = zeros(1,size(Ypp,2));
    
    for i=1:kfold
        %     (round(size(Ypp,1)/kfold))
        tst = (1+(i-1):kfold:size(Ypp,1));
        Xtstc = Xpp;
        Xtstc(tst,:) = [];
        Xtstp = Xpp(tst,:);
        Ytstc = Ypp;
        Ytstc(tst,:) = [];
        for j = 1:size(Ytstc,2)
            [~,~,~,~,b(:,j),~,~,~] = plsregress(Xtstc,Ytstc(:,j),ncomp(j));
            Ytstp(:,j) = [ones(size(Xtstp,1),1),Xtstp]*b(:,j);
            Ycv(tst,j) = Ytstp(:,j);
            Ytstp = [];
        end
    end
    for j = 1:size(Ycal,2)
        if ppy == 1
            Ycv(:,j) = Ycv(:,j) + mean(Ycal(:,j));
        elseif ppy == 2
            Ycv(:,j) = Ycv(:,j).*std(Ycal(:,j)) + mean(Ycal(:,j));
        elseif ppy == 3
            Ycv(:,j) = Ycv(:,j).*sqrt(std(Ycal(:,j))) + mean(Ycal(:,j));
        end
        
        rmse(2,j) = sqrt(mean((Ycal(:,j)-Ycv(:,j)).^2));
        r = (corrcoef(Ycal(:,j),Ycv(:,j))).^2;
        Rcv(1,j) = r(1,2);
        R(2,:) = Rcv;
    end    
end

    
%% Residuals and Hotteling
waitbar(7/9,h1,['PLS calculation: ', sprintf('%0.1u',7) '/', sprintf('%0.1u',9)]);

%Residuals
if ptype == 1 % PLS1
    
    for i=1:size(stats_final.Xresiduals,1)
        Q(i,1) = stats_final.Xresiduals(i,:)*stats_final.Xresiduals(i,:)';
    end
    
    %Qcrit
    [m,n] = size(stats_final.Xresiduals);
    if m>n
        [~,s,~] = svd(stats_final.Xresiduals);s = diag(s);
    else
        [~,s,~] = svd(stats_final.Xresiduals');s = diag(s);
    end
    s = s.^2/(m-1);
    % s = Q;
    cl = 0.95;
    pc = 0;
    theta1 = sum(s(pc+1:end,1));
    theta2 = sum(s(pc+1:end,1).^2);
    theta3 = sum(s(pc+1:end,1).^3);
    h0 = 1-2*theta1*theta3/3/(theta2.^2);
    ca    = sqrt(2)*erfinv(2*cl-1);
    h0    = ca*sqrt(2*theta2*h0.^2)/theta1;
    h2    = theta2*h0*(h0-1)/(theta1.^2);
    Qcrit = theta1*(1+h0+h2).^(1/h0);
    
    Qind = find(Q>Qcrit);
    
    %Hotteling
    Fcrit = finv(0.95,size(Xpp,1),ncomp);
    Tcrit = (ncomp*((size(Xpp,1)^2)-1)/(size(Xpp,1)*(size(Xpp,1)-ncomp))*Fcrit);
    t2ind = find(stats_final.T2>Tcrit);
    
else % PLS2
    
    for j=1:size(Ycal,2)
        for i=1:size(stats_final{j}.Xresiduals,1)
            Q(i,j) = stats_final{j}.Xresiduals(i,:)*stats_final{j}.Xresiduals(i,:)';
        end
        
        %Qcrit
        [m,n] = size(stats_final{j}.Xresiduals);
        if m>n
            [~,s,~] = svd(stats_final{j}.Xresiduals);s = diag(s);
        else
            [~,s,~] = svd(stats_final{j}.Xresiduals');s = diag(s);
        end
        s = s.^2/(m-1);
        cl = 0.95;
        pc = 0;
        theta1 = sum(s(pc+1:end,1));
        theta2 = sum(s(pc+1:end,1).^2);
        theta3 = sum(s(pc+1:end,1).^3);
        h0 = 1-2*theta1*theta3/3/(theta2.^2);
        ca    = sqrt(2)*erfinv(2*cl-1);
        h0    = ca*sqrt(2*theta2*h0.^2)/theta1;
        h2    = theta2*h0*(h0-1)/(theta1.^2);
        Qcrit(:,j) = theta1*(1+h0+h2).^(1/h0);
        
        Qind{j} = find(Q(:,j)>Qcrit(:,j));
        
        %Hotteling
        Fcrit = finv(0.95,size(Xpp,1),ncomp(j));
        Tcrit(j) = (ncomp(j)*((size(Xpp,1)^2)-1)/(size(Xpp,1)*(size(Xpp,1)-ncomp(j)))*Fcrit);
        t2ind{j} = find(stats_final{j}.T2>Tcrit(j));
        
    end
end


%% Report
waitbar(7/9,h1,['PLS calculation: ', sprintf('%0.1u',7) '/', sprintf('%0.1u',9)]);
% 
% Report_cal(:,1) = {'Sample number', 'Latent Variables', 'Variables number', 'RMSEC', 'RMSECV', 'R²', 'R²CV', 'Minimum Y', 'Maximum Y', '% Error'};
% Report_cal(:,1) = {'Sample number', 'Latent Variables', 'Variables number', 'RMSEC', 'RMSECV', 'R²', 'R²CV', 'Minimum Y', 'Maximum Y', '% Error', 'Bench version'}; TO BE IMPLMENTED
% Report_cal{1,2} = size(Xpp,1);
% Report_cal{2,2} = ncomp;
% Report_cal{3,2} = size(Xpp,2);
% Report_cal{4,2} = rmse(1);
% Report_cal{5,2} = rmse(2);
% Report_cal{6,2} = R(1,2);
% Report_cal{7,2} = Rcv(1,2);
% Report_cal{8,2} = min(Yc);
% Report_cal{9,2} = max(Yc);
% Report_cal{10,2} = rmse(2)/(max(Yc) - min(Yc))*100;
% Report_cal{11,2} = bench version TO BE IMPLEMENTED
% 
% xlswrite('Report_cal.xlsx', Report_cal);

Model.beta = BETA;
Model.Xloading = XL;
Model.Xscores = XS;
Model.LV = ncomp;
Model.ppx = ppx;
Model.ppy = ppy;
Model.stats = stats_final;
if ptype == 1 % PLS1
    Model.stats.rmse = rmse;
    Model.stats.R(1,:) = R;
    Model.stats.R(2,:) = Rcv;
    Model.stats.Tcrit = Tcrit;
    Model.stats.Qcrit = Qcrit;
    Model.stats.Q = Q;
    Model.stats.T2 = stats_final.T2;
else
    for i=1:size(Ycal,2)
        Model.stats{i}.rmse = rmse(:,i);
        Model.stats{i}.R(1,:) = R(1,i);
        Model.stats{i}.R(2,:) = Rcv(1,i);
        Model.stats{i}.Tcrit = Tcrit(i);
        Model.stats{i}.Qcrit = Qcrit(i);
        Model.stats{i}.Q = Q(:,i);
        Model.stats{i}.T2 = stats_final{i}.T2;
    end
end
Model.Yfit = Yc;
Model.Ycv = Ycv;
Model.Yc = Ycal;
Model.mX = mean(Xcal);
Model.sX = std(Xcal - ones(size(Xcal,1),1) * mean(Xcal));
% Model.Xcal_new = Xcal(:,:);
% save ('calibration_model.mat','Model');

%% Permutation
% Not implemented yet
waitbar(8/9,h1,['PLS calculation...: ', sprintf('%0.1u',8) '/', sprintf('%0.1u',9)]);
% 
% if exist('options','var')
%     if isfield(options,'nPerms')
%         [R2_perm,~] = plsper_LQA(Xpp,Ypp,ncomp,options.nPerms,options);
%     end
% end

%% Plots
waitbar(9/9,h1,['PLS calculation...: ', sprintf('%0.1u',9) '/', sprintf('%0.1u',9)]);

if plots >= 1
    for y = 1:1:size(Ycal,2)
        %Calibration figure
        axes1 = axes('Parent',figure);
        sc1=scatter(Ycal(:,y),Yc(:,y),40,'filled'); hold on;
        axesLimits1 = xlim(axes1);
        ax = axis;
        xplot1 = linspace(ax(1),ax(2));
        fitResults1 = polyfit(Ycal(:,y),Yc(:,y), 1);
        yplot1 = polyval(fitResults1, xplot1); hold on;
        fitLine1 = plot(xplot1,yplot1,'DisplayName','   linear','Parent',axes1,...
            'Tag','linear',...
            'Color',[0.47 0.67 0.19],'linestyle','--');
        %     Rfl = refline(1,0);
        %     Rfl.Color = [0.91 0.41 0.10];
        %     Rfl.LineStyle = '-';
        ax = axis;
        [~,im] = max(ax);
        Rfl = line([ax(im-1) ax(im)],[ax(im-1) ax(im)],'Color',[0.91 0.41 0.10],'linestyle','-');
        %     legend('Calibration Samples','Adjusted Line','Reference Line')
        xlabel 'Reference Y'
        ylabel 'Predicted Y'
        title 'Calibration correlation'
        annotation('textbox',...
            [0.16 0.73 0.3 0.15],...
            'String',{['RMSEC = ' num2str(rmse(1,y))],['R² = ' num2str(R(1,y))]},...
            'FitBoxToText','off',...
            'EdgeColor','none');
        box 'on'
        grid 'on'
        
        %CV figure
        axes2 = axes('Parent',figure);
        sc2=scatter(Ycal(:,y),Ycv(:,y),40,'filled'); hold on;
        axesLimits2 = xlim(axes2);
        ax2 = axis;
        xplot2 = linspace(ax2(1),ax2(2));
        fitResults2 = polyfit(Ycal(:,y),Ycv(:,y), 1);
        yplot2 = polyval(fitResults2, xplot2); hold on;
        fitLine2 = plot(xplot2,yplot2,'DisplayName','   linear','Parent',axes2,...
            'Tag','linear',...
            'Color',[0.47 0.67 0.19],'linestyle','--');
        ax2 = axis;
        [~,im] = max(ax2);
        Rfl = line([ax2(im-1) ax2(im)],[ax2(im-1) ax2(im)],'Color',[0.91 0.41 0.10],'linestyle','-');
        %     Rfl = refline(1,0);
        %     Rfl.Color = [0.91 0.41 0.10];
        %     Rfl.LineStyle = '-';
        %     legend('Calibration Samples', 'Adjusted Line','Reference Line')
        xlabel 'Reference Y'
        ylabel 'Predicted CV Y'
        title 'CV correlation'
        annotation('textbox',...
            [0.16 0.73 0.3 0.15],...
            'String',{['RMSECV = ' num2str(rmse(2,y))],['R²CV = ' num2str(Rcv(1,y))]},...
            'FitBoxToText','off',...
            'EdgeColor','none');
        box 'on'
        grid 'on'
    end
end

if ptype == 1
    q = Qind;
    Qind = cell(1,1);
    Qind{1} = q;
    
    st_final = stats_final;
    stats_final = cell(1,1);
    stats_final{1} = st_final;

    t2i = t2ind;
    t2ind = cell(1,1);
    t2ind{1} = t2i;
end
if plots == 2
    for y=1:size(Ycal,2)
        %DModX figure
        axes1 = axes('Parent',figure);
        plot(Q(:,y),'-o','MarkerFaceColor','auto'); hold on;
        if ~isempty(Qind{y})
            scatter(Qind{y},Q(Qind{y},y),40,'filled');hold on
            text(Qind{y},Q(Qind{y},y),Am(Qind{y}), 'Interpreter', 'none');
        end
        ax = axis;
        line([ax(1) ax(2)],[Qcrit(y) Qcrit(y)],'Color','k','LineStyle','--');axis tight
        ylim ([ax(3)*1.2 ax(4)*1.2])
        xlabel 'Samples'
        ylabel 'Distance to the model in X'
        title 'Residuals'
        box 'on'        
        
        %Hotteling figure
        axes1 = axes('Parent',figure);
        plot(stats_final{y}.T2,'-or','MarkerFaceColor','auto'); hold on;
        if ~isempty(t2ind{y})
            scatter(t2ind{y},stats_final{y}.T2(t2ind{y}),40,'r','filled');hold on
            text(t2ind{y},stats_final{y}.T2(t2ind{y}),Am(t2ind{y}), 'Interpreter', 'none');
        end
        ax = axis;
        line([ax(1) ax(2)],[Tcrit(y) Tcrit(y)],'Color','k','LineStyle','--'); axis tight;
        ylim ([ax(3)*1.2 ax(4)*1.2])
        xlabel 'Samples'
        ylabel 'T^2'
        title 'Hotteling T^2'
        box 'on'
    end
end

close (h1)
