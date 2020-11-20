function [p_treat,p_treat_norm,rmspe_treat_log,rmspe_treat_diff,diff_treat_log_norm,corr_treat] = calc_p_treat(series_mat,series_cand,mc_cand,treat_time,v_opt,timevec,c_treat,hypothesis,excl_outlier,b_plot)

% CALC_P_TREAT Calculates placebo treatments and p-values of significant 
% differences for the treated.
% Authors: Makram El-Shagi, Axel Lindner and Gregor von Schweinitz
% _________________________________________________________________________
%
% SYNTAX:
%
% [p_treat,p_treat_norm,rmspe_treat_log,rmspe_treat_diff,diff_treat_log_norm,corr_treat] = calc_p_treat(series_mat,series_cand,mc_cand,treat_time,v_opt,timevec,c_treat,hypothesis,excl_outlier,b_plot)
% _________________________________________________________________________
%
% INPUT
% series_mat            matrix of series (true and synthetic)
% series_cand           variable of interest in candidate countries
% mc_cand               matching criteria in candidate countries
% treat_time            timing of the treatment in the series
% v_opt                 optimal v from synthetic matching
% timevec               timing of observations
% c_treat               treatment countries 
% hypothesis            vector of hypothesis for one-sided tests (see below)
% excl_outlier          boolean if outliers should be excluded
% b_plot                boolean if p-values should be plotted
% _________________________________________________________________________
%
% OUTPUT
%
% p_treat               p-value of over-/undervaluation per period, not
%                       normalized
% p_treat_norm          p-value of over-/undervaluation per period,
%                       normalized by pre-treatment standard deviation
% rmspe_treat_log       rmspe of treatments (log) 
% rmspe_treat_diff      rmspe of treatments (abs) 
% diff_treat_log_norm   errors of log-observable, normalized by rmspe
% corr_treat            correlation of treatment with observable using
%                       split samples

if nargin<10
    b_plot = 0;
end
if nargin<9
    excl_outlier = 0;
end
if nargin<8            % general hypothesis: 1 - country is overvalued; -1 - country is undervalued; 0 - period by period one-sided test
    hypothesis = zeros(size(series_mat,2),1);
end


%% init
[T,n_treat] = size(series_mat);
st = series_mat(:,1:2:n_treat-1);
sts = series_mat(:,2:2:n_treat);
dt_log = log(st)-log(sts);
dtp = (st-sts)./st;
dtpn = dtp;

n_cand = size(series_cand,2);
if excl_outlier
    n_cand=n_cand-2;
end
n_treat = n_treat/2;
p_treat = zeros(T,n_treat);
p_treat_norm = zeros(T,n_treat);
rmspe_treat_log = zeros(1,n_treat);
rmspe_treat_diff = zeros(1,n_treat);
corr_treat = zeros(2,n_treat);
[times,~,tind] = unique(treat_time);
mc_use = mc_cand(1:length(v_opt),:);

fac = (n_cand-1)/(n_cand);
add = 1/(2*n_cand);

%% calc period-by-period tests
for k = 1:length(times)
    post = find(tind==k);
    std_t = nanstd(dtp(1:times(k),post));
    
    rmspe_treat_diff(post) = std_t;
    rmspe_treat_log(post) = nanstd(dt_log(1:times(k),post));
    corr_treat(1,post) = diag(corr(st(1:times(k),post),sts(1:times(k),post),'rows','pairwise'));
    corr_treat(2,post) = diag(corr(st(times(k)+1:end,post),sts(times(k)+1:end,post),'rows','pairwise'));
    
    dtpn(:,post) = dtp(:,post)./(repmat(std_t,T,1));
    [~,~,dcp,dcpn] = calc_candtreat(series_cand,mc_use,times(k),v_opt);
    if excl_outlier
        dcp_mean = mean(dcp);
        pmax = find(dcp_mean==max(dcp_mean));
        pmin = find(dcp_mean==min(dcp_mean));
        dcp(:,[pmin pmax]) = [];
        dcpn(:,[pmin pmax]) = [];
    end

    for c = 1:length(post)
        switch hypothesis(post(c))
            case 0
                overval = dtp(:,post(c))>0;
                p_treat(overval,post(c)) = 1-nanmean((repmat(dtp(overval,post(c)),1,n_cand)-dcp(overval,:))>0,2);
                p_treat(~overval,post(c)) = 1-nanmean((repmat(dtp(~overval,post(c)),1,n_cand)-dcp(~overval,:))<0,2);

                p_treat_norm(overval,post(c)) = 1-nanmean((repmat(dtpn(overval,post(c)),1,n_cand)-dcpn(overval,:))>0,2);
                p_treat_norm(~overval,post(c)) = 1-nanmean((repmat(dtpn(~overval,post(c)),1,n_cand)-dcpn(~overval,:))<0,2);

            case 1
                p_treat(:,post(c)) = 1-nanmean((repmat(dtp(:,post(c)),1,n_cand)-dcp)>0,2);
                p_treat_norm(:,post(c)) = 1-nanmean((repmat(dtpn(:,post(c)),1,n_cand)-dcpn)>0,2);
            case -1
                p_treat(:,post(c)) = 1-nanmean((repmat(dtp(:,post(c)),1,n_cand)-dcp)<0,2);
                p_treat_norm(:,post(c)) = 1-nanmean((repmat(dtpn(:,post(c)),1,n_cand)-dcpn)<0,2);
        end
        
        %adjustment of nan
        posnan = find(isnan(dtp(:,post(c))));
        for p = 1:length(posnan)
            p_treat(posnan(p),post(c)) = p_treat(posnan(p)-1,post(c));
            p_treat_norm(posnan(p),post(c)) = p_treat_norm(posnan(p)-1,post(c));
        end
        
        %adjustment of extreme values
    end
end
diff_treat_log_norm = dt_log./repmat(rmspe_treat_log,T,1);

p_treat = fac*p_treat+add;
p_treat_norm = fac*p_treat_norm+add;

%% plot period-by-period tests
if b_plot
    scrsz = get(0,'ScreenSize');
    scrsz(2) = 31;
    scrsz(4) = scrsz(4)-30;
    h = figure('OuterPosition',scrsz);
    set(h,'defaulttextinterpreter','latex');
    ncols = ceil(sqrt(n_treat));
    nrows = ceil(n_treat/ncols);
    pleg = [0.05 0.01 0.9 0.05];

    b1 = 0;
    b2 = 1;

    t = timevec;
    if iscell(timevec)
        t_plot = datevec(timevec,'dd.mm.yyyy');
        t = datenum(t_plot);
        t = t-t_plot(:,3)+1;            %bringing it to the first of the month
        years = (t_plot(1,1):t_plot(end,1));
        t_ticks = nan(length(years),1);
        for k = 1:length(years)
            t_ticks(k) = find(t_plot(:,1)==years(k),1);
        end

        %October 2008
        t_cr = find((t_plot(:,1)==2008) .* (t_plot(:,2)==10));
    else
        t_cr = find(t==2008);
    end


    for k= 1:n_treat

        h = subplot(nrows,ncols,k);
        fill([t(t_cr) t(t_cr:end)' t(end)],[b1 repmat(b2,1,length(t(t_cr:end))) b1],[0.8 0.8 0.8],'EdgeColor','none');
        hold on
        plot([t(treat_time(k)) t(treat_time(k))],[b1 b2],'k','LineWidth',2);
        plot(t,[p_treat(:,k) p_treat_norm(:,k)]);
        plot(t,0.2*ones(length(t),1),'k')
        if iscell(timevec)
            set(gca,'XTick',t_ticks)
            datetick('x','yy')
        end
        axis tight

        xlabel('Year','FontSize',16)
        ylabel('Probability','FontSize',16)
        title(c_treat(k),'FontSize',20)
        if k == n_treat      
            lh=legend(h,{'Crisis period';'Introduction of the Euro';'Probability of being significantly over/undervalued';'Normed Probability'},'FontSize',12,'Location','NorthWest');
            set(lh,'Orientation','horizontal')
            set(lh,'OuterPosition',pleg);
        end
    end
end