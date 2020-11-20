function [ssr_tot,ssr,w_mat] = ev_v_quadprog_panel(v,series_treat,series_cand,mc_treat,mc_cand,treat_time_vec,opt_qp)

% EV_V_QUADPROG_PANEL evaluates a v for the synthetic matching.
% Called by synthesize_fmincon_panel
% Source: Abadie/Gardezabal (2003)
% Authors: Makram El-Shagi, Axel Lindner and Gregor von Schweinitz
% _________________________________________________________________________
%
% SYNTAX:
%
% [ssr_tot,ssr,w_mat] = ev_v_quadprog_panel(v,series_treat,series_cand,mc_treat,mc_cand,treat_time_vec,opt_qp)
% _________________________________________________________________________
%
% INPUT
% 
% v                 importance vector v to be evaluated
% series_treat      is an t x n_treat matrix with the treated series
% series_cand       is an t x n_cand matrix with the series of candidates
% mc_treat          is an m x n_treat matrix with the criteria to be met
%                   (treated)
% mc_cand           is an m x n_cand matrix with the criteria to be met
%                   (candidates)
% treat_time_vec    is the timing of the treatment in the series
% opt_qp            option set for quadratic programming
% _________________________________________________________________________
%
% OUTPUT
%
% ssr_tot           total sum of squared residuals (over all treated)
% ssr               sum of squared residuals (per treated individuum)
% w_mat             n_cand x n_treat matrix of weights for candidates

if nargin<7
    opt_qp = optimset('quadprog');
end

n_treat = size(series_treat,2);
n_cand = size(series_cand,2);
w_mat = zeros(n_cand,n_treat);
ssr = zeros(n_treat,1);

X0 = mc_cand;
D = diag(v);
H = X0'*D*X0;
H = (H+H')/2;
for n = 1:n_treat
    X1 = mc_treat(:,n);
    f = - X1'*D*X0;
    [w]=quadprog(H,f,[],[],ones(1,n_cand),1,zeros(n_cand,1),ones(n_cand,1),[],opt_qp);
    w(w<0)=0;
    w = w/sum(w);
    w_mat(:,n) = w;
    
    Z1 = series_treat(1:treat_time_vec(n)-1,n);    
    Z0 = series_cand(1:treat_time_vec(n)-1,:);
    one_zero=isnan(Z0);
    w_rep=repmat(w',size(Z0,1),1);
    w_rep(one_zero)=0;
    weight_adjustment=sum(w_rep,2).^(-1);
    weight=w_rep .* repmat(weight_adjustment,1,size(Z0,2));
    Z0(one_zero)=0;
    synthetic = sum(Z0.*weight,2);

    e = Z1 - synthetic;
    %should we accept a bias? For function minimization, we probably should...
    e = e-nanmean(e);

    ssr(n) = nansum(e.^2);
end

ssr_tot = sum(ssr);
end

