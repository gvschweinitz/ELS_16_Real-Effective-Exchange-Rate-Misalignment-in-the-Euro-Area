function [mc_treat,mc_cand,mc_est,mc_coeff,mc_latent,nanpos,pc_options] = prepare_estimation_panel(mc_treat_orig,mc_cand_orig,pc,normalization)

% PREPARE_ESTIMATION_PANEL constructs the matching-matrix and estimation
% series for Abadie-Gardeazabal (2003)
% Authors: Makram El-Shagi, Axel Lindner and Gregor von Schweinitz
% _________________________________________________________________________
%
% SYNTAX:
%
% [mc_treat,mc_cand,mc_est,mc_coeff,mc_latent,nanpos,pc_options] = 
%           prepare_estimation(mc_treat_orig,mc_cand_orig,pc,normalization)
% _________________________________________________________________________
%
% INPUT
%
% mc_treat_orig:        M x N_1 matrix with the matching criteria of treatment countries
% mc_cand_orig:         M x N_0 matrix with the matching criteria of candidate countries
% pc:                   logical variable if matching is to be performed on
%                           principal components of matching variables
% normalization:        logical variable if data are to be normalized
% _________________________________________________________________________
%
% OUTPUT
%
% mc_treat:             M x N_1 matrix of matching criteria (or pc) of treatment
%                           countries after adjustment
% mc_cand:              M x N_0 matrix of matching criteria (or pc) of treatment
%                           countries after adjustment
% mc_est:               M x (N_1+N_0) matrix of all matching criteria (or pc)
% mc_coeff:             M x M matrix of loadings: mc_coeff * mc = mc_pc
% mc_latent:            M x 1 vector of eigenvalues
% nanpos:               positions of NaN in original matching criteria
% pc_options:           Structure with mc_latent, and selection criterion
%                           for number of components
if nargin<4
    normalization = 1;
end
if nargin<3
    pc = 0;
end
if (~normalization) && pc
    warning('Data have to be normalized for principal components');
    normalization = 1;
end

n_treat = size(mc_treat_orig,2);
mc = [mc_treat_orig mc_cand_orig];
nanpos = sum(isnan(mc_treat_orig),2)>0;
mc(nanpos==1,:)=[];
if normalization
    %standard-normalization
    mc = mc';
    std(mc);
    mc_m = nanmean(mc);
    mc_std = nanstd(mc);
    mc_rows = size(mc,1);
    mc = (mc-repmat(mc_m,mc_rows,1))./repmat(mc_std,mc_rows,1);
    mc_est = mc;
    mc_coeff = eye(size(mc,2));
    mc_latent = ones(size(mc,2),1);

    %principal components
    if pc == 1
        [mc_coeff, mc_est, mc_latent] = pca(mc);
    end

    mc_est = mc_est';
    mc_coeff = mc_coeff';
else
    mc_est = mc;
    mc_coeff = eye(size(mc,1));
    mc_latent = ones(size(mc,1));
end

mc_treat = mc_est(:,1:n_treat);
mc_cand = mc_est(:,n_treat+1:end);
pc_options = struct('mc_latent',mc_latent,'type','ev','limit',1);