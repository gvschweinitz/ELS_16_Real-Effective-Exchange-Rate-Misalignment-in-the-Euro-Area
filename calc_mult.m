function calc_mult(series_treat,series_cand,mc_treat,mc_cand,treat_time_vec,runs,startrun,filename,parproc,pc,pc_options)

% CALC_MULT calculates the optimal synthetic matching from multiple random
% starting values. Random starting values are needed in cases where the
% optimization problem of synthetic matching has a non-standard surface.
% Runs are divided into blocks (based on the number of parallel
% processors), and the results of each block are saved separately. This
% allows for combining the results of many consecutive calls of calc_mult,
% for example in case of limited computing power.
% Authors: Makram El-Shagi, Axel Lindner and Gregor von Schweinitz
% _________________________________________________________________________
%
% SYNTAX:
% calc_mult(series_treat,series_cand,mc_treat,mc_cand,treat_time_vec,runs,startrun,filename,parproc,pc,pc_options)
% _________________________________________________________________________
%
% INPUT
% series_treat      TxN_1 matrix of variable of interest for treatment countries
% series_cand       TxN_0 matrix of variable of interest for candidate countries
% mc_treat          MxN_1 matrix of matching criteria for treatment countries
% mc_cand           MxN_0 matrix of matching criteria for candidate countries
% treat_time_vec    time position (index) of treatment time
% runs              Number of runs with random starting values
% startrun          First run to be done (can be larger than one)
% filename          .mat-file where results are to be saved
% parproc           Number of parallel processors
% pc                Boolean if principal components of matching criteria should be used
% pc_options:       Structure with principal components of matching criteria, and selection criterion
%                       for number of components
% 
% _________________________________________________________________________
%
% OUTPUT (SAVED IN FILENAME)
% w_runs:           Country weights, by run
% v_runs:           Weights of matching criteriy, by run
% ssr_runs:         Value of optimization function, by run

if nargin<10
    pc = 0;
end
if nargin<11
    if pc
        error('need to have structure of pc options');
    end
end

blockmult = 40;
if pc
    fprintf('Limiting matching criteria... ')
    if strcmp(pc_options.type, 'ev')
        mc_use = pc_options.mc_latent>pc_options.limit;
        fprintf('using %4.0f of %4.0f\n',sum(mc_use),size(mc_treat,1));
    elseif strcmp(pc_options.type, 'perc')
        share = cumsum(pc_options.mc_latent)/sum(pc_options.mc_latent);
        mc_use = [1:find(share>pc_options.limit,1)];
        fprintf('using %4.0f of %4.0f\n',max(mc_use),size(mc_treat,1));
    end
    
    mc_treat = mc_treat(mc_use,:);
    mc_cand = mc_cand(mc_use,:);
    size(mc_treat)
    size(mc_cand)
end

if startrun>1
    out = open([filename '.mat']);
    w_runs = out.w_runs;
    v_runs = out.v_runs;
    ssr_runs = out.ssr_runs;
    clear out
else
    w_runs = nan(runs,size(mc_cand,2),size(mc_treat,2));
    v_runs = nan(runs,size(mc_cand,1));
    ssr_runs = nan(runs,size(mc_treat,2));
end

n_outer = ceil((runs-startrun)/(parproc*blockmult));

for j = 1:n_outer

    k1 = min(startrun + (j-1)*parproc*blockmult,runs);
    k2 = min(startrun + j*parproc*blockmult - 1,runs);
    fprintf('Starting block %4.0f\n',j)
    tic
    parfor (k = k1:k2,parproc)
        c=clock;
        rand(ceil(c(6)),1); %randomizing starting value parfor
        fprintf('\t\t\t starting run %4.0f\n',k)
        [w_mat, v, ssr]=synthesize_fmincon_panel(series_treat,series_cand,mc_treat,mc_cand,treat_time_vec);
        w_runs(k,:,:) = w_mat;
        v_runs(k,:) = v;
        ssr_runs(k,:) = ssr;
    end
    ssr_min = min(sum(ssr_runs(1:k2,:),2));
    fprintf('\t ---- block %4.0f finished (runs %4.0f to %4.0f; minimum ssr until now %4.2f; time taken %4.0f sec) ---- \n',j,k1,k2,ssr_min,toc);
    save(filename,'w_runs','v_runs','ssr_runs','-append')    
end
