function  plot_and_save_bw(series_treat,series_cand,timevec,w_mat,treat_time,c_treat,plotcontrol,savenme,plotpos)

% PLOT_AND_SAVE_BW plots and saves Figures 3 and 4 of 
% Real Effective Exchange Rate Misalignment in the Euro Area: A
% Counterfactual Analysis
% Authors: Makram El-Shagi, Axel Lindner and Gregor von Schweinitz
% _________________________________________________________________________
%
% SYNTAX:
% plot_and_save_bw(series_treat,series_cand,timevec,w_mat,treat_time,c_treat,plotcontrol,savenme,plotpos)
% _________________________________________________________________________
%
% INPUT
% series_treat      TxN_1 matrix of variable of interest for treatment countries
% series_cand       TxN_0 matrix of variable of interest for candidate countries
% timevec           Tx1 string cell of observation times
% w_mat             N_0xN_1 matrix of optimal country weights
% treat_time        time position (index) of treatment time
% c_treat           names of treatment countries
% plotcontrol       potential subselection of countries to be plotted
% savenme           file name of the figure to be saved
% plotpos           position of subplot of every treatment country

b_save = 1;
if nargin<8
    b_save=0;
end

if nargin <7
    nc = length(c_treat);
    plotcontrol = (1:nc);
else
    nc = length(plotcontrol);
end

if nargin<9
    plotpos = 1:nc;
end

b1 = 60;
b2 = 140;
t = timevec;
if iscell(timevec)
    t_plot = datevec(timevec,'mm/dd/yyyy');
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

series_treat = series_treat(:,plotcontrol);
c_treat = c_treat(plotcontrol);
treat_time = treat_time(plotcontrol);
w_mat = w_mat(:,plotcontrol);
T = length(timevec);
synth = zeros(T,nc);
Z1 = zeros(T,nc);
one_zero=isnan(series_cand);
series_cand(one_zero) = 0;
for k=1:nc
    w = w_mat(:,k)';
    w_repZ=repmat(w,T,1);
    w_repZ(one_zero)=0;
    w_adjZ=sum(w_repZ,2).^(-1);
    weightZ=w_repZ .* repmat(w_adjZ,1,size(series_cand,2));
    synth(:,k) = sum(series_cand .* weightZ,2);
    e = series_treat(1:treat_time(k)-1,k)-synth(1:treat_time(k)-1,k);
    synth(:,k) = exp(synth(:,k) + nanmean(e));
    
    denom = nanmean(exp(series_treat(1:treat_time(k)-1,k)));
    Z1(:,k) = exp(series_treat(:,k))/denom*100;
    synth(:,k) = synth(:,k)/denom*100;
end


scrsz = get(0,'ScreenSize');
scrsz(2) = 31;
scrsz(4) = scrsz(4)-30;
h = figure('OuterPosition',scrsz);
set(h,'defaulttextinterpreter','latex');

ncols = ceil(sqrt(nc));
nrows = ceil(nc/ncols);

for k = 1:nc
    h = subplot(nrows,ncols,plotpos(k));
    pleg = [0.05 0.01 0.9 0.05];
    fill([t(t_cr) t(t_cr:end)' t(end)],[b1 repmat(b2,1,length(t(t_cr:end))) b1],[0.8 0.8 0.8],'EdgeColor','none');
    hold on
    plot([t(treat_time(k)) t(treat_time(k))],[b1 b2],'k','LineWidth',2);
    line_fewer_markers(t,Z1(:,k),length(t_ticks),'k-o','LineWidth',1.2,'MarkerSize',6);
    line_fewer_markers(t,synth(:,k),length(t_ticks),'k--*','LineWidth',1.2,'MarkerSize',6);
    if iscell(timevec)
        set(gca,'XTick',t_ticks)
        datetick('x','yy')
    end
    axis tight
    title(c_treat(k),'FontSize',20)
    if k == nc      
        lh=legend(h,{'Crisis period';'Introduction of the Euro';'observed';'synthetic'},'FontSize',16,'Location','NorthWest');
        set(lh,'Orientation','horizontal')
        set(lh,'Position',pleg);
    end
end

if b_save
    saveas2([savenme '.pdf'],300,'pdf')
end