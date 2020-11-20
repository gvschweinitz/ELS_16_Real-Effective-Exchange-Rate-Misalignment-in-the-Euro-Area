function [series_cand_treat,diff_cand_log_norm,diff_cand_perc,diff_cand_perc_norm] = calc_candtreat(series_cand,mc_cand,treat_time,v_opt)

% CALC_CANDTREAT calculates placebo synthetic matching for candidate
% countries using optimal weights V for matching criteria. In the
% estimation, candidate countries are taken in turn as placebo treatment
% country, with all other candidates forming the group for the synthetic
% (placebo) treatment).
% Authors: Makram El-Shagi, Axel Lindner and Gregor von Schweinitz
% _________________________________________________________________________
%
% SYNTAX:
% [series_cand_treat,diff_cand_log_norm,diff_cand_perc,diff_cand_perc_norm] = calc_candtreat(series_cand,mc_cand,treat_time,v_opt)
% _________________________________________________________________________
%
% INPUT
% series_cand       TxN_0 matrix of variable of interest for candidate
%                       countries (in logs)
% mc_cand           MxN_0 matrix of matching criteria for candidate countries
% treat_time        time position (index) of treatment time. Only the first
%                       element is used in case this is a vector
% 
% _________________________________________________________________________
%
% OUTPUT
% series_cand_treat:    Variable of interest of countries under placebo treatment
% diff_cand_log_norm:   Difference between observed and synthetic series in
%                           logs, divided by pre-treatment estimation error
% diff_cand_perc:       Difference between observed and synthetic series in
%                           levels, divided by observed series (in levels)
% diff_cand_perc_norm:  Difference between observed and synthetic series in
%                           levels, divided by observed series and estimation errors (in levels)

if size(mc_cand,2)>length(v_opt)
    mc_cand = mc_cand(1:length(v_opt),:);
end
[T,n_cand] = size(series_cand);
series_cand_treat = zeros(T,n_cand);

treat_time = treat_time(1);

for k = 1:n_cand
    sct = series_cand(:,k);
    scc = series_cand;
    scc(:,k)=[];
    
    mct = mc_cand(:,k);
    mcc = mc_cand;
    mcc(:,k)=[];
    
    [~,~,w] = ev_v_quadprog_panel(v_opt,sct,scc,mct,mcc,treat_time);
    
    one_zero=isnan(scc);
    w_rep=repmat(w',size(scc,1),1);
    w_rep(one_zero)=0;
    weight_adjustment=sum(w_rep,2).^(-1);
    weight=w_rep .* repmat(weight_adjustment,1,size(scc,2));
    scc(one_zero)=0;
    series_cand_treat(:,k) = sum(scc.*weight,2);

end

e = series_cand(1:treat_time-1,:) - series_cand_treat(1:treat_time-1,:);
series_cand_treat = series_cand_treat + repmat(nanmean(e),T,1);
diff_cand_log_norm = (series_cand-series_cand_treat)./repmat(std(e),T,1);
diff_cand = (exp(series_cand)-exp(series_cand_treat));
diff_cand_perc = diff_cand./exp(series_cand);
e = diff_cand_perc(1:treat_time-1,:);
diff_cand_perc_norm = diff_cand_perc./repmat(std(e),T,1);
