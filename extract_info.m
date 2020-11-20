function [series_treat,series_cand,mc_treat_orig,mc_cand_orig,time,treat_time,choice,c_treat,c_cand,countries_long,indicators_criteria,eurocountries,matchingcountries,controls,mc_avail] = ...
    extract_info(file,treat_date,controls,b_log)

% EXTRACT_INFO extracts all necessary data for running a (panel version of)
% synthetic matching from an Excel file. If the file does not have the same
% format as "Data_REERmonthly_level_032014.xlsx", this code may need
% adaptation.
% Authors: Makram El-Shagi, Axel Lindner and Gregor von Schweinitz
% _________________________________________________________________________
%
% SYNTAX:
%
% [series_treat,series_cand,mc_treat_orig,mc_cand_orig,time,treat_time,choice,c_treat,c_cand,countries_long,indicators_criteria,eurocountries,matchingcountries,controls,mc_avail] = ...
%    extract_info(file,treat_date,controls,b_log)
% _________________________________________________________________________
%
% INPUT
%
% file:         name of the excel file (incl ".xlsx") containing all necessary data.
%                   The first two sheets contain control group assignments 
%                   of countries and summary statistics to be derived from
%                   matching criteria. Summary statistics can be mean, std,
%                   growth rates and five-year averages
%                   The third sheet contains the variable of interest.
%                   remaining sheets contain original data on matching criteria
% treat_date:   N_1 x 1 vector with treatment years (treatment in January of the year)
% controls:     Structure of options for construction of matching criteria. Includes:
%                   min_mean: min #obs for calculation of mean (default: 1)
%                   min_std: min #obs for calculation of std (default: 2)
%                   min_growth: min #obs for calculation of (av) growth
%                   rate (default: 5)
%                   min_5year: min #obs for 5-year averages (default: 1)
%                   max_countries: max #candidate countries where a
%                   matching criterion is allowed to be missing. If more,
%                   matching criterion is dropped
% b_log:        Boolean if variable of interest needs log-transformation
% _________________________________________________________________________
%
% OUTPUT
%
% series_treat:     TxN_1 matrix of variable of interest in treatment countries
% series_cand:      TxN_0 matrix of variable of interest in candidate countries
% mc_treat_orig:    M x N_1 matrix with the matching criteria of treatment countries
% mc_cand_orig:     M x N_0 matrix with the matching criteria of candidate countries
% time:             time vector of series_treat/series_cand    
% treat_time:       N_1x1 vector of indices of treatment times in "time"
% choice:               positions of NaN in original matching criteria
% c_treat:          ISO-2 names of treatment countries
% c_cand:           ISO-2 names of candidate countries
% countries_long:   Names of all countries in the data file
% indicators_criteria: Names of all matching criteria
% eurocountries:    Boolean vector, 1 if treatment country
% matchingcountries:Boolean vector, 1 if candidate country
% controls:         Structure of final controls (incl defaults)
% mc_avail:         Number of observations used for the calculation of the
%                       matching criterion.

% file should contain '.xls' or 'xlsx' depending on the type of file.
if nargin<5
    b_log = 1;
end
if nargin<4
    controls = struct();         
end
if ~isfield(controls,'min_mean')
    min_mean = 1;
    controls.min_mean = min_mean;
else
    min_mean = controls.min_mean;
end
if ~isfield(controls,'min_std')
    min_std = 2;
    controls.min_std = min_std;
else
    min_std = controls.min_std;
end
if ~isfield(controls,'min_growth')
    min_growth = 5;
    controls.min_growth = min_growth;
else
    min_growth = controls.min_growth;
end
if ~isfield(controls,'min_5year')
    min_5year = 1;
    controls.min_5year= min_5year;
else
    min_5year = controls.min_5year;
end
if ~isfield(controls,'max_countries')
    % maximum number of countries missing in matchingcountries and
    % eurocountries for the criterion to be used
    max_countries = 10;
    controls.max_countries= max_countries;
else
    max_countries = controls.max_countries;
end

[~, sheets] = xlsfinfo(file);
indicators = sheets(3:end);
n_indicators = length(indicators);

% ControlSheet_Country
[data,text] = xlsread(file,sheets{1});
countries = text(1,2:end);
countries_long = text(2,2:end);
eurocountries = data(4,:);
matchingcountries = data(2,:);
n_countries = length(countries);

% ControlSheet_Ind
[~,text] = xlsread(file,sheets{2});
ind_check = text(1,2:end);
if length(ind_check)~=length(indicators)
    error('wrong indicators')
end
for k = 1:length(indicators)
    if ~strcmp(indicators{k},ind_check{k})
        error('wrong indicators')
    end
end

choice=text(2:end,2:end);
s_prepost = zeros(size(choice,1),1);
s_transformation = zeros(size(choice,1),1);
for k = 1:size(choice,1)
    comp = lower(text{k+1,1});
    if isequal(comp,'prepost')
        s_prepost(k) = 1;
    elseif isequal(comp,'transformation')
        s_transformation(k) = 1;
    end
end
choice_time = choice(s_prepost==1,:);
choice_transformations = choice(s_transformation==1,:);
el_times = size(choice_time,1);
el_transforms=size(choice_transformations,1);

% Initialization of Matching Criteria
count=0;
for k = 1:n_indicators
    count_time = 0;
    count_trans = 0;
    for k_time = 1:size(choice_time,1)
        if ~isempty(choice_time{k_time,k})
            count_time = count_time+1;
        end
    end
    for k_trans = 1:size(choice_transformations,1)
        if ~isempty(choice_transformations{k_trans,k})
            count_trans = count_trans+1;
        end
    end
    count = count + count_time*count_trans;
end
matching_criteria = nan(count,n_countries);
indicators_criteria = cell(count,1);
mc_avail = nan(count,n_countries);

% loading the matching data
[data,text] = xlsread(file,sheets{3});
for c = 1:n_countries
    if ~strcmp(countries{c},text{1,c+1})
        error('wrong country')
    end
end
t = size(data,1);
date_min  = [];
if size(data,2)==n_countries+1
    % Date column in data matrix
    if isnumeric(data(1,1))
        % Yearly data
        t_min = data(1,1);
    else
        % Monthly data
        date_min = datevec(data(1,1),'mm/dd/yyyy');
        time = data(:,1);
    end
    series = data(:,2:end);
elseif size(data,2)==n_countries
    % Date column in text cell
    date_min = datevec(text{2,1},'mm/dd/yyyy');
    time = text(2:t+1,1);
    series = data;
else
    error('Very bad sizing of REER')
end

first_event = min(treat_date);
last_event = max(treat_date);
treat_time = zeros(length(treat_date),1);
if ~isempty(date_min)
    t_min = date_min(1);
    treat_mc_first = first_event-t_min+1;
    treat_mc_last = last_event-t_min+1;
    dvec = datevec(time,'mm/dd/yyyy');
    for k = 1:length(treat_date)
        treat_time(k) = find(treat_date(k) == dvec(:,1),1);
    end
else
    treat_mc_first = first_event-t_min+1;
    treat_mc_last = last_event-t_min+1;
    treat_time = treat_date-t_min+1;
end

counter = 1;
for n = 2:n_indicators
    [data_ind,text] = xlsread(file,sheets{2+n});
    for c = 1:n_countries
        if ~strcmp(countries{c},text{1,c+1})
            countries{c}
            indicators{n}
            error('wrong country')
        end
    end
    if size(data_ind,2) == n_countries+1
        data = data_ind(:,2:end);
        timevec = data_ind(:,1);
    else
        data = data_ind;
        timevec = str2num(cell2mat(text(2:end,1)));
    end
%     if size(data,1)~=t
%         error('not enough years')
%     end
    for m_time = 1:el_times
        if ~isempty(choice_time{m_time,n})
            option_time=lower(choice_time{m_time,n});
            for m_trans = 1:el_transforms
                if ~isempty(choice_transformations{m_trans,n})
                    option_trans=lower(choice_transformations{m_trans,n});
                    fprintf('%s ---- %s --- %s\n', indicators{n},option_time,option_trans)
                    
                    x = nan(1,n_countries);
                    if isequal(option_time,'pre')
                        data_select=data(1:treat_mc_first-1,:);
                        time_select = timevec(1:treat_mc_first-1);
                    end
                    if isequal(option_time,'post')
                        data_select=data(treat_mc_last:end,:);
                        time_select = timevec(treat_mc_last:end);
                    end
                    if size(data_select,1) ~= size(time_select)
                        size(data_select,1)
                        size(time_select)
                    end
                    
                    if isequal(option_trans,'mean')
                        cols = sum(~isnan(data_select))>=min_mean;
                        x(cols) = nanmean(data_select(:,cols));
                    end
                    if isequal(option_trans,'median')
                        cols = sum(~isnan(data_select))>=min_median;
                        x(cols) = nanmean(data_select(:,cols));
                    end
                    if isequal(option_trans,'std')
                        cols = sum(~isnan(data_select))>=min_std;
                        if n==2 %GDP
                            growth = diff(log(data_select(:,cols)));
                            x(cols) = nanstd(growth);
                        else
                            x(cols) = nanstd(data_select(:,cols));
                        end
                    end
                    if isequal(option_trans,'growth')
                        for c = 1:n_countries
                            notnan = find(~isnan(data_select(:,c)));
                            if isempty(notnan)
                                countries{c}
                                indicators{n}
                            else

                                first = notnan(1);
                                last = notnan(end);

                                if (last-first) >= min_growth
                                    x(c) = (log(data_select(last,c)) - log(data_select(first,c))) / (last-first);
                                end
                            end
                        end
                    end
                    if isequal(option_trans,'5year')
                        T = size(data_select,1);
                        for b = 1:ceil(T/5)
                            x = nan(1,n_countries);
                            if isequal(option_time,'pre')
                                pos = max(T-b*5+1,1):T-(b-1)*5;
                            else
                                pos = (b-1)*5+1:min(b*5,T);
                            end
                            block = [num2str(time_select(pos(1))) '-' num2str(time_select(pos(end)))];
                            cols = sum(~isnan(data_select(pos,:)))>=min_5year;
                            if sum(cols)>0
                                x(cols) = nanmean(data_select(pos,cols));
                            end
                            
                            matching_criteria(counter,:)=x;
                            indicators_criteria{counter} = [indicators{n} '-' option_time '_' option_trans '_' block];
                            mc_avail(counter,:) = sum(~isnan(data_select(pos,:)));
                            counter=counter+1;
                        end
                    end
                    

%                         if isnan(x)
%                             fprintf('Criterion %s/%s in country %s not available\n',option_time,option_trans,countries{c})
%                         end
                    if ~isequal(option_trans,'5year')
                        matching_criteria(counter,:)=x;
                        indicators_criteria{counter} = [indicators{n} '-' option_time '_' option_trans];
                        mc_avail(counter,:) = sum(~isnan(data_select));
                        counter=counter+1;
                    end
                end
            end
        end
    end
end
mc_missing_eur = sum(mc_avail(:,eurocountries == 1)==0,2);
mc_missing_tot = sum(mc_avail(:,eurocountries == 1 | matchingcountries==1)==0,2);
pos_use = (mc_missing_eur==0) & (mc_missing_tot<=max_countries);
fprintf('Deleting criteria \n')
indicators_criteria(~pos_use)
matching_criteria = matching_criteria(pos_use,:);
indicators_criteria = indicators_criteria(pos_use);
mc_avail = mc_avail(pos_use,:);

series_treat = series(:,eurocountries==1);
series_cand = series(:,matchingcountries==1);
mc_treat_orig = matching_criteria(:,eurocountries==1);
mc_cand_orig = matching_criteria(:,matchingcountries==1);
c_treat = countries(eurocountries==1);
c_cand = countries(matchingcountries==1);
datafile_name = file;

if b_log
    series_treat = log(series_treat);
    series_cand = log(series_cand);
end

end


