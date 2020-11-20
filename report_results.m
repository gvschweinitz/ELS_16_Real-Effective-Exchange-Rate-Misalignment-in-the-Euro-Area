function [series_mat,avdev,mc_mat,descr_stat,descr_corr,v_pc,v_mc,w_mat] = report_results(series_treat,series_cand,mc_treat_orig,mc_cand_orig,mc_coeff,ssr_runs,w_runs,v_runs,treat_time,crisis_time,nanpos,ordering)

% REPORT_RESULTS combines data with estimation output to identify the
% optimal solution of the synthetic matching, and derive the output of
% Tables 2-5 of 
% Real Effective Exchange Rate Misalignment in the Euro Area: A
% Counterfactual Analysis
% Authors: Makram El-Shagi, Axel Lindner and Gregor von Schweinitz
% _________________________________________________________________________
%
% SYNTAX:
% [series_mat,avdev,mc_mat,descr_stat,descr_corr,v_pc,v_mc,w_mat] = ...
%       report_results(series_treat,series_cand,mc_treat_orig,mc_cand_orig,mc_coeff,ssr_runs,w_runs,v_runs,treat_time,crisis_time,nanpos,ordering)
% _________________________________________________________________________
%
% INPUT
% series_treat      TxN_1 matrix of variable of interest for treatment countries
% series_cand       TxN_0 matrix of variable of interest for candidate countries
% mc_treat_orig:    M x N_1 matrix with the matching criteria of treatment countries
% mc_cand_orig:     M x N_0 matrix with the matching criteria of candidate countries
% mc_coeff:             M x M matrix of loadings: mc_coeff * mc = mc_pc
% ssr_runs:         Value of optimization function, by run
% w_runs:           Country weights, by run
% v_runs:           Weights of matching criteriy, by run
% treat_time        N_1x1 vector (index) of treatment time
% crisis_time       Index: start of the financial crisis
% nanpos:               positions of NaN in original matching criteria
% ordering          N_0x1 vector ordering of control countries
% 
% _________________________________________________________________________
%
% OUTPUT
% series_mat        Tx(2*N_1) matrix of series (true and synthetic)
% avdev:            5xN_1 matrix of average differences (Table 4, plus a
%                   row for the deviation at the start of the financial
%                   crisis)
% mc_mat:           Mx(2*N_1) matrix of matching criteria (true and
%                   synthetic, Table 6)
% descr_stat:       Descriptive statistics of matching criteria (Table 2)
% descr_corr:       Correlation of matching criteria (Table 3)
% v_pc:             V-Weights on principal components (Table 5)
% v_mc:             MxM matrix V-Weights on matching criteria. Not needed
% w_mat:            N_0xN_1 matrix of optimal country weights

nanord = [];
if nargin<12
    ordering = (1:size(mc_treat_orig,1));
else
    if length(ordering) > size(mc_treat_orig,1)
        if sum(~isnan(ordering)) ~= size(mc_treat_orig,1)
            error('wrong ordering')
        else
            nind = length(ordering);
            nanord = isnan(ordering);
            ordering = ordering(~nanord);
        end
    end
end

pos = find(sum(ssr_runs,2)==min(sum(ssr_runs,2)),1);
w_mat = squeeze(w_runs(pos,:,:));

[T,nc] = size(series_treat);
synth = zeros(T,nc);
one_zero = isnan(series_cand);
series_cand(one_zero) = 0;
for k=1:nc
    w = w_mat(:,k)';
    w_repZ = repmat(w,T,1);
    w_repZ(one_zero) = 0;
    w_adjZ = sum(w_repZ,2).^(-1);
    weightZ = w_repZ .* repmat(w_adjZ,1,size(series_cand,2));
    var = sum(series_cand .* weightZ,2);
    e = series_treat(:,k) - var;
    var = var + nanmean(e(1:treat_time(k)));
    synth(:,k) = var;
end


series_mat = zeros(T,nc);
mc_mat = zeros(size(mc_treat_orig));

synth = exp(synth);
Z1 = exp(series_treat);
mc_synth = mc_cand_orig*w_mat;
e = (Z1-synth)./Z1;
avdev = zeros(5,nc);
for k = 1:nc
    avdev(:,k) = [nanmean(e(1:treat_time(k)-1,k)); nanmean(e(treat_time(k):crisis_time-1,k)); nanmean(e(crisis_time:end,k)); e(crisis_time-1,k); nanmean(e(crisis_time-11:crisis_time,k))];
    series_mat(:,(k-1)*2+1) = Z1(:,k);
    series_mat(:,k*2) = synth(:,k);
    mc_mat(:,(k-1)*2+1) = mc_treat_orig(:,k);
    mc_mat(:,k*2) = mc_synth(:,k);
end

descr_stat = [nanmean(mc_treat_orig,2) nanstd(mc_treat_orig,[],2) nanmin(mc_treat_orig,[],2) nanmax(mc_treat_orig,[],2) nanmin(mc_cand_orig,[],2) nanmax(mc_cand_orig,[],2)];
descr_corr = corr([mc_treat_orig mc_cand_orig]','rows','pairwise');

v = v_runs(pos,:);
v(v<0) = 0;
v_mc = nan(size(mc_treat_orig,1),size(mc_treat_orig,1));
v_mc(nanpos==0,nanpos==0) = mc_coeff(1:numel(v),:)'*diag(v)*mc_coeff(1:numel(v),:);
v_pc = v';

mc_mat = mc_mat(ordering,:);
descr_stat = descr_stat(ordering,:);
descr_corr = descr_corr(ordering,ordering);
v_mc = v_mc(ordering,ordering);

if ~isempty(nanord)
    temp = nan(nind,size(mc_mat,2));
    temp(nanord==0,:) = mc_mat;
    mc_mat = temp;
    
    temp = nan(nind,size(descr_stat,2));
    temp(nanord==0,:) = descr_stat;
    descr_stat = temp;
    
    temp = nan(nind,nind);
    temp(nanord==0,nanord==0) = descr_corr;
    descr_corr = temp;

    temp = nan(nind,nind);
    temp(nanord==0,nanord==0) = v_mc;
    v_mc = temp;
end