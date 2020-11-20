function plot_placebo_bw(diff_cand,diff_treat,timevec,treat_time,savenme)

% PLOT_PLACEBO_BW plots and saves Figure 2 of
% Real Effective Exchange Rate Misalignment in the Euro Area: A
% Counterfactual Analysis
% Authors: Makram El-Shagi, Axel Lindner and Gregor von Schweinitz
% _________________________________________________________________________
%
% SYNTAX:
% plot_placebo_bw(diff_cand,diff_treat,timevec,treat_time,savenme)
% _________________________________________________________________________
%
% INPUT
% diff_cand         TxN_0 matrix; Difference between observed and synthetic
%                       series, candidate countries
% diff_treat        TxN_1 matrix; Difference between observed and synthetic
%                       series, treatment countries
% timevec           Tx1 string cell of observation times
% treat_time        time position (index) of treatment time
% savenme           file name of the figure to be saved

b_save=1;
if nargin<5
    b_save = 0;
end

scrsz = get(0,'ScreenSize');
scrsz(2) = scrsz(2)+30;
scrsz(4) = scrsz(4)-30;
h = figure('OuterPosition',scrsz);

treat_time = unique(treat_time);

b1 = -6;
b2 = 6;

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

set(gca,'FontSize',12)
subplot(2,1,1);
hdl = 0;
hdl(1) = fill([t(t_cr) t(t_cr:end)' t(end)],[b1 repmat(b2,1,length(t(t_cr:end))) b1],[0.8 0.8 0.8],'EdgeColor','none');
hold on
for k = 1:length(treat_time)
    hdl(k+1) = plot([t(treat_time(k)) t(treat_time(k))],[b1 b2],'k','LineWidth',2);
end
plot(t,diff_treat,'k')
if iscell(timevec)
    set(gca,'XTick',t_ticks)
    datetick('x','yy')
end
title('(a) Normalized prediction errors, treatment countries (Euro introduction in 1999 and 2001)','FontSize',12)
axis tight
xlabel('Year','FontSize',12)
ylabel('Prediction error / RMSPE(1980-treatment)','FontSize',12)

subplot(2,1,2);
hdl = 0;
hdl(1) = fill([t(t_cr) t(t_cr:end)' t(end)],[b1 repmat(b2,1,length(t(t_cr:end))) b1],[0.8 0.8 0.8],'EdgeColor','none');
hold on
plot([t(treat_time(1)) t(treat_time(1))],[b1 b2],'k','LineWidth',2);
plot(t,diff_cand,'k')
if iscell(timevec)
    set(gca,'XTick',t_ticks)
    datetick('x','yy')
end
title('(b) Normalized prediction errors, placebo studies (Euro introduction in 1999)','FontSize',12)
legend({'Crisis period';'Introd. of the Euro'},'Location','southoutside','Orientation','horizontal');
axis tight
xlabel('Year','FontSize',12)
ylabel('Prediction error / RMSPE(1980-1999)','FontSize',12)

hold off
if b_save
    saveas2([savenme '.pdf'],300,'pdf')
end