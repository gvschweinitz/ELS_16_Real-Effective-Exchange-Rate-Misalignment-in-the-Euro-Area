function [w_mat, v, ssr]=synthesize_fmincon_panel(series_treat,series_cand,mc_treat,mc_cand,treat_time_vec)

% SYNTHESIZE_FMINCON_PANEL applies the synthetic matching estimator. This 
% is the panel-version of the program.
% Calls ev_v_quadprog_panel
% Source: Abadie/Gardezabal (2003)
% Authors: Makram El-Shagi, Axel Lindner and Gregor von Schweinitz
% _________________________________________________________________________
%
% SYNTAX:
%
% [weights, V, synthetic] 
%           = synthesize(series,mc,synth_country,treat_time)
% _________________________________________________________________________
%
% INPUT
%
% series_treat      is an t x n_treat matrix with the treated series
% series_cand       is an t x n_cand matrix with the series of candidates
% mc_treat          is an m x n_treat matrix with the criteria to be met
%                   (treated)
% mc_cand           is an m x n_cand matrix with the criteria to be met
%                   (candidates)
% treat_time_vec    is the timing of the treatment in the series
% _________________________________________________________________________
%
% OUTPUT
%
% w_mat             n_cand x n_treat matrix of weights for candidates
% v                 contains the diagonal entries of the (diag) importance
%                      matrix (m x m) of criteria
% ssr               sum of squared residuals for series


opt_fmin = optimset('fmincon');
opt_fmin = optimset(opt_fmin,'MaxFunEvals',10000,'display','none');

opt_qp = optimset('quadprog');
opt_qp = optimset(opt_qp,'display','none');

warning off


%sum of v normalized to unity
l = size(mc_treat,1);
v_start = rand(1,l);
v_start = v_start/sum(v_start);
% v_start = 1/l * ones(1,l);
% fprintf('Starting optimization')
[v_opt] = fmincon(@(v) ev_v_quadprog_panel(v,series_treat,series_cand,mc_treat,mc_cand,treat_time_vec,opt_qp),...
    v_start,[],[],ones(1,l),1,zeros(l,1),[],[],opt_fmin);

% fprintf('Optimization finished')

[~,ssr,w_mat] = ev_v_quadprog_panel(v_opt,series_treat,series_cand,mc_treat,mc_cand,treat_time_vec,opt_qp);
if size(v_opt,1) == l-1
    v=[1 v_opt];
else
    v = v_opt;
end

warning on

