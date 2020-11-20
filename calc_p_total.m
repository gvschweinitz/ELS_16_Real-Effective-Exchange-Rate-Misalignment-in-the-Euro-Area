function [p_smooth,p_tot] = calc_p_total(p_mat,treat_time,crisis_time,window,control_test)
% Aggregates p-values from many tests to deal with multiple testing problem
% use calc_p_Bonferoni, calc_p_Zstat or calc_p_Fisher.
% Authors: Makram El-Shagi, Axel Lindner and Gregor von Schweinitz
% _________________________________________________________________________
% INPUT
% p_mat:        individual p-values, dimension TxN
% treat_time:   Treatment time (start testing there)
% crisis_time:  Time of crisis as a separator for larger blocks
% window:       Window length for smoothing
% control_test: String to switch between tests
% _________________________________________________________________________
% OUTPUT
% p_smooth:     aggregated p-values over blocks of periods
% p_tot:        aggregated p-values over three blocks, plus two information rows: 
%                   (1) full evaluation period
%                   (2) Pre-crisis 
%                   (3) Crisis
%                   (4) Input p-values one period before the crisis
%                   (5) Smoothed p-values one period before the crisis

[T,N] = size(p_mat);
p_smooth = nan(T,N);
p_tot = nan(5,N);

if nargin<5
    control_test = 'Bonferoni';
end
alpha_test = [0.01 0.05 0.1];

for k = 1:N
    for t = treat_time(k)+window-1:T
        switch control_test
            case 'Fisher'
                p_smooth(t,k) = calc_p_Fisher(p_mat(t-window+1:t,k));
            case 'Z-stat'
                p_smooth(t,k) = calc_p_Zstat(p_mat(t-window+1:t,k));
            case 'Bonferoni'
                p_smooth(t,k) = calc_p_Bonferoni(p_mat(t-window+1:t,k),alpha_test);
        end
    end
    switch control_test
        case 'Fisher'
            p_tot(1,k) = calc_p_Fisher(p_mat(1:treat_time(k)-1,k));
            p_tot(2,k) = calc_p_Fisher(p_mat(treat_time(k):crisis_time-1,k));
            p_tot(3,k) = calc_p_Fisher(p_mat(crisis_time:end,k));
        case 'Z-stat'
            p_tot(1,k) = calc_p_Zstat(p_mat(1:treat_time(k)-1,k));
            p_tot(2,k) = calc_p_Zstat(p_mat(treat_time(k):crisis_time-1,k));
            p_tot(3,k) = calc_p_Zstat(p_mat(crisis_time:end,k));
        case 'Bonferoni'
            p_tot(1,k) = calc_p_Bonferoni(p_mat(1:treat_time(k)-1,k),alpha_test);
            p_tot(2,k) = calc_p_Bonferoni(p_mat(treat_time(k):crisis_time-1,k),alpha_test);
            p_tot(3,k) = calc_p_Bonferoni(p_mat(crisis_time:end,k),alpha_test);
    end
    
    p_tot(4,k) = p_mat(crisis_time-1,k);
    p_tot(5,k) = p_smooth(crisis_time-1,k);
end