function plot_pie(X,str_label,str_title,col_select,ordering,savenme)

% PLOT_PIE plots and saves Figures 6 and 7 of 
% Real Effective Exchange Rate Misalignment in the Euro Area: A
% Counterfactual Analysis
% Authors: Makram El-Shagi, Axel Lindner and Gregor von Schweinitz
% _________________________________________________________________________
%
% SYNTAX:
% plot_pie(X,str_label,str_title,col_select,ordering,savenme)
% _________________________________________________________________________
%
% INPUT
% X                 N_0xN_1 matrix of weights, sum by column to one.
% str_label         N_0x1 vector of string labels for rows in X (pie parts)
% str_title         N_1x1 vector of labels for columns in X (subplots)
% col_select        potential subselection of countries to be plotted
% ordering          Ordering of candidate countries (for color selection)
% savenme           file name of the figure to be saved

b_save = 1;
if nargin<6
    b_save=0;
end
if nargin<5
    ordering = 1:size(X,1);
end

[~,sortpos] = sort(ordering);
X = X(sortpos,:);
str_label = str_label(sortpos);
ordering = ordering(sortpos);

cutoff = 0.05;
norm = 1;

nplot = length(col_select);
ncols = ceil(sqrt(nplot));
nrows = ceil(nplot/ncols);

if size(str_label,1)<size(str_label,2)
    str_label = str_label';
end
scrsz = get(0,'ScreenSize');
scrsz(2) = 31;
scrsz(4) = scrsz(4)-30;
h = figure('OuterPosition',scrsz);


pos_cutoff = sum(X>cutoff,2);

if min(ordering)==-1 && max(ordering) == 1
    c = zeros(length(ordering),3);
    c(ordering==1,:) = repmat([0 0.8 0],sum(ordering==1),1);
    c(ordering==-1,:) = repmat([1 0 0],sum(ordering==-1),1);
    c(ordering==0,:) = repmat([1 1 0.1],sum(ordering==0),1);
    c = [c(pos_cutoff>0,:);[0 0.8 0]];
else
    c = colormap(jet);
    ncolor = sum(pos_cutoff>0);
    cpos = 1:(size(c,1)-1)/(ncolor):size(c,1);
    cpos = round(cpos);
    c = c(cpos,:);
    c=c(ordering(1:ncolor+1),:);
end
Xplot = X(pos_cutoff>0,col_select);
size(str_label)
size(pos_cutoff)
str_label = str_label(pos_cutoff>0);
tit_plot = str_title(col_select);
for k = 1:nplot
    hs = subplot(nrows,ncols,k);
    Xpos = (Xplot(:,k)>=cutoff);
    Xtemp = Xplot(Xpos,k);
    Xtemp = [Xtemp;norm-sum(Xtemp)];
    lab_temp = [str_label(Xpos);'Other'];
%     lab_temp = strcat(lab_temp,': ');
    c_temp = [c(Xpos,:);c(end,:)];
    
    hp = pie(Xtemp);
    set(hp(2:2:end),'FontSize',20)
    
    %work on texts
    hText = findobj(hp,'Type','text'); % text handles
    percentValues = get(hText,'String'); % percent values
    
    if length(percentValues)<length(lab_temp)
        lab_temp=lab_temp(1:end-1);
        c_temp = c_temp(1:end-1,:);
    end
    combinedstrings = strcat(lab_temp,': ',percentValues); % text and percent values
    
    oldExtents_cell = get(hText,'Extent'); % cell array
    oldExtents = cell2mat(oldExtents_cell); % numeric array
    set(hText,{'String'},combinedstrings);
    newExtents_cell = get(hText,'Extent'); % cell array
    newExtents = cell2mat(newExtents_cell); % numeric array
    width_change = newExtents(:,3)-oldExtents(:,3);
    signValues = sign(oldExtents(:,1));
    offset = signValues.*(1.1*width_change/2);
    textPositions_cell = get(hText,{'Position'}); % cell array
    textPositions = cell2mat(textPositions_cell); % numeric array
    textPositions(:,1) = textPositions(:,1) + offset; % add offset
    set(hText,{'Position'},num2cell(textPositions,2)) % set new position
    
    colormap(c_temp);
    title(tit_plot(k),'FontSize',24);
    
    freezeColors
end

if b_save
    saveas2([savenme '.pdf'],300,'pdf')
end