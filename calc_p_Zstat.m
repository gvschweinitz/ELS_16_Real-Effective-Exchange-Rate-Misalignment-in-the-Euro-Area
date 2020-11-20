function p = calc_p_Zstat(p_vec)

% CALC_P_ZSTAT calculates the Z-statistic for a series of p-values. If
% p_vec is a matrix, it calculates the Z-statistic by column.
% Authors: Makram El-Shagi, Axel Lindner and Gregor von Schweinitz
% _________________________________________________________________________
%
% SYNTAX:
% p = calc_p_Zstat(p_vec)
% _________________________________________________________________________
%
% INPUT
% p_vec             Vector of observation-specific p-values
% _________________________________________________________________________
%
% OUTPUT
% p                 Z-statistic (by column)

p_vec(p_vec==0) = eps;
p_vec(p_vec==1) = 1-eps;

test = sum(norminv(p_vec,0,1))/sqrt(length(p_vec));
p = normcdf(test,0,1);