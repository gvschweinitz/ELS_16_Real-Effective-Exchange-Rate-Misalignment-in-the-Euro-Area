function p = calc_p_Fisher(p_vec)

% CALC_P_FISHER runs a Fisher-test for joint significance.
% Attention: the Fisher test assumes independence of p-values
% Source: Maddala-Wu (1999)
% Authors: Makram El-Shagi, Axel Lindner and Gregor von Schweinitz
% _________________________________________________________________________
%
% SYNTAX:
% p = calc_p_Fisher(p_vec)
% _________________________________________________________________________
%
% INPUT
% p_vec             Vector of observation-specific p-values
% _________________________________________________________________________
%
% OUTPUT
% p                 joint p-value


p_vec(p_vec==0) = eps;
p_vec(p_vec==1) = 1-eps;

test = -2*sum(log(p_vec));
p = 1-chi2cdf(test,2*length(p_vec));