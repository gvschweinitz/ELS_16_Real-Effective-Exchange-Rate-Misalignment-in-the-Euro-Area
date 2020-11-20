% Main code to replicate baseline results for the paper
% Real Effective Exchange Rate Misalignment in the Euro Area: A
% Counterfactual Analysis
% by Makram El-Shagi, Axel Lindner and Gregor von Schweinitz 
% (Review of International Economics, 24(1), 37-66, 2016)

clear;

%% Parameters
% Files
datafile_name = 'Data_REERmonthly_level_032014.xlsx';
filename = 'calc_Benchmark.mat';
% Treatment variable specification
treat_date = [1999;1999;1999;1999;1999;2001;1999;1999;1999;1999;1999];
b_log = 1;                  % Boolean if REER should be in logs
% Country classification
c_core = [1:5 9];           % Indices of core countries
c_per = [6:8 10 11];        % Indices of periphery countries
ordering_pie = [-1;1;-1;-1;-1;1;0;-1;1;1;-1;-1;-1;1;1;0;-1;1;1;-1;1;1;-1]; %ordering of candidate countries (developed/developing)
% Control variables
controls = struct('min_mean',1,'min_std',5,'min_growth',5);
pc = 1;                     % Boolean if principal components of matching variables should be used
% Calculation settings
parproc = 16;                % Number of parallel processors to be used for estimation
runs = 20000;                 % Number of random starting values of V
% Settings for tests
excl_outliers = 0;          % Boolean if two most extreme placebo countries should be excluded in calculation of placebo treatments
p_hypo = [-1;-1;-1;-1;-1;1;1;1;-1;1;1];
save(filename)

%% Core code
% Prepare Calculations
[series_treat,series_cand,mc_treat_orig,mc_cand_orig,time,treat_time,choice,c_treat,c_cand,countries_long,indicators_criteria,eurocountries,matchingcountries,controls,mc_avail] = extract_info(datafile_name,treat_date,controls,b_log);
[mc_treat,mc_cand,mc_est,mc_coeff,mc_latent,nanpos,pc_options] = prepare_estimation_panel(mc_treat_orig,mc_cand_orig,pc,1);
ordering = 1:size(mc_treat_orig,1);
ordering = ordering(nanpos==0);
crisis_time = find(strcmp(time,'10/15/2008'));

% Model estimation
save(filename)
calc_mult(series_treat,series_cand,mc_treat,mc_cand,treat_time,runs,1,filename,parproc,pc,pc_options);  % Output saved in filename
load(filename)

% Report results
[series_mat,avdev,mc_mat,descr_stat,descr_corr,v_pc,v_mc,w_mat] = report_results(series_treat,series_cand,mc_treat_orig,mc_cand_orig,mc_coeff,ssr_runs,w_runs,v_runs,treat_time,crisis_time,nanpos,ordering);
[series_cand_treat,diff_cand_log_norm,diff_cand_perc,diff_cand_perc_norm] = calc_candtreat(series_cand,mc_cand,treat_time,v_pc);
[p_treat,p_treat_norm,rmspe_treat_log,rmspe_treat_diff,diff_treat_log_norm,corr_treat] = calc_p_treat(series_mat,series_cand,mc_cand,treat_time,v_pc,time,c_treat,p_hypo,excl_outliers,0);
[p_smooth_norm_Bonf,p_tot_norm_Bonf] = calc_p_total(p_treat_norm,treat_time,crisis_time,12,'Bonferoni');
[p_smooth_norm_Fisher,p_tot_norm_Fisher] = calc_p_total(p_treat_norm,treat_time,crisis_time,12,'Fisher');
save(filename)



%% Plots and tables
load('calc_Benchmark.mat')
c_treatl = countries_long(eurocountries==1);
% Plots in black and white. Adapt function for color version
plot_REER_dev_bw(series_treat,'longrun',time,treat_time,c_treatl,b_log,[5 4 9 6 7 11],'Fig1_bw')
plot_placebo_bw(diff_cand_log_norm,diff_treat_log_norm,time,treat_time,'Fig2_bw')                               
plot_and_save_bw(series_treat,series_cand,time,w_mat,treat_time,c_treatl,c_core,'Fig3_bw')
plot_and_save_bw(series_treat,series_cand,time,w_mat,treat_time,c_treatl,c_per,'Fig4_bw')

% Data for Figure 5 -> call function plotavdev.m
[~,sorting] = sort(avdev(5,:));
temp = cell(11,3);
temp(:,1) = c_treatl(sorting);
temp(:,2) = num2cell(avdev(5,sorting))';
temp(1:5,3) = {'Undervaluation / Fair Valuation'};
temp(6:8,3) = {'Modest Overvaluation'};
temp(9:11,3) = {'Overvaluation'};
temp = cell2table(temp,'VariableNames',{'Country';'Deviation';'Label'});
writetable(temp,'plotavdev.csv')  

% Data for Table 7
temp = array2table([w_mat ordering_pie],'VariableNames',[c_treatl {'Class'}],'RowNames',c_cand); 
writetable(temp,'w_mat.csv','WriteRowNames',true)

plot_pie(w_mat,c_cand,c_treatl,c_core,ordering_pie,'wpie_core_noTrade') % Figure 6
plot_pie(w_mat,c_cand,c_treatl,c_per,ordering_pie,'wpie_per_noTrade')   % Figure 7

%% Robustness checks. Need to be calculated separately
% clear all
% load('calc_robust_5year.mat')
% plot_and_save_bw(series_treat,series_cand,time,w_mat,treat_time,c_treatl,1:11,'Fig8_bw')
% 
% clear all
% load('calc_robust_noTrade.mat')
% plot_and_save_bw(series_treat,series_cand,time,w_mat,treat_time,c_treatl,1:11,'Fig9_bw')
% 
% clear all
% load('calc_robust_Inflation.mat')
% plot_and_save_bw(series_treat,series_cand,time,w_mat,treat_time,c_treatl,1:11,'Fig10_bw')
% 
% clear all
% load('calc_robust_noBG.mat')
% plot_and_save_bw(series_treat,series_cand,time,w_mat,treat_time,c_treatl,1:10,'Fig11_bw',[1 3:11])
