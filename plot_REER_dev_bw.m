function plot_REER_dev_bw(series_treat,norm,timevec,treat_time,c_treatl,b_log,c_sel,savenme)

% PLOT_REER_DEV_BW plots and saves Figure 1 of
% Real Effective Exchange Rate Misalignment in the Euro Area: A
% Counterfactual Analysis
% Authors: Makram El-Shagi, Axel Lindner and Gregor von Schweinitz
% _________________________________________________________________________
%
% SYNTAX:
% plot_REER_dev_bw(series_treat,norm,timevec,treat_time,c_treatl,b_log,c_sel,p_hypo,savenme)
% _________________________________________________________________________
%
% INPUT
% series_treat      TxN_1 matrix of variable of interest in treatment countries
% norm              'Euro' normalize series to be 100 at treatment date
%                   'longrun' normalize series to be 100 on average
%                   pre-treatment
% timevec           Tx1 string cell of observation times
% treat_time        time position (index) of treatment time
% c_treatl          names of treatment countries
% b_log             Boolean if variable of interest needs log-transformation
% c_sel             potential subselection of countries to be plotted        
% savenme           file name of the figure to be saved

[T,N] = size(series_treat);
b_save = 1;
if nargin<8
    b_save = 0;
end
if nargin<7
    c_sel = 1:N;
end

series = series_treat;
if b_log
    series = exp(series);
end

denom = ones(N,1);
switch norm
    case 'Euro'
        for n = 1:N
            denom(n) = series(treat_time(n),n);
        end
    case 'longrun'
        for n = 1:N
            denom(n) = nanmean(series(1:treat_time(n)-1,n));
        end
end

for n=1:N
    series(:,n) = series(:,n)/denom(n)*100;
end
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

b1 = 60;
b2 = 140;
series = series(:,c_sel);
c_treatl = c_treatl(c_sel);
treat_time = unique(treat_time);

scrsz = get(0,'ScreenSize');
scrsz(2) = 31;
scrsz(4) = scrsz(4)-30;
h = figure('OuterPosition',scrsz);
set(h,'defaulttextinterpreter','latex');

hdl = 0;
hdl(1) = fill([t(t_cr) t(t_cr:end)' t(end)],[b1 repmat(b2,1,length(t(t_cr:end))) b1],[0.8 0.8 0.8],'EdgeColor','none');
hold on
hdl(2) = line_fewer_markers(t,series(:,1),length(t_ticks),'k-o','LineWidth',1.3,'MarkerSize',8);
hdl(3) = line_fewer_markers(t,series(:,2),length(t_ticks),'k-.o','LineWidth',1.3,'MarkerSize',8);
hdl(4) = line_fewer_markers(t,series(:,3),length(t_ticks),'k--o','LineWidth',1.3,'MarkerSize',8);
hdl(5) = line_fewer_markers(t,series(:,4),length(t_ticks),'k-*','LineWidth',1.3,'MarkerSize',8);
hdl(6) = line_fewer_markers(t,series(:,5),length(t_ticks),'k-.*','LineWidth',1.3,'MarkerSize',8);
hdl(7) = line_fewer_markers(t,series(:,6),length(t_ticks),'k--*','LineWidth',1.3,'MarkerSize',8);

for k = 1:length(treat_time)
    hdl(k+length(c_sel)+1) = plot([t(treat_time(k)) t(treat_time(k))],[b1 b2],'k','LineWidth',2);
end

if iscell(timevec)
    set(gca,'XTick',t_ticks)
    datetick('x','yy')
end
xlabel('Time','FontSize',16);
ylabel('Real Effective Exchange Rate','FontSize',16);
legend(['Crisis period' c_treatl 'Introduction of the Euro'],'FontSize',16,'Location','NorthWest')
hold off
if b_save
    saveas2([savenme '.pdf'],300,'pdf')
end