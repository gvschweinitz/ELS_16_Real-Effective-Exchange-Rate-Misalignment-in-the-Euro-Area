function p = calc_p_Bonferoni(p_vec,alpha)

% CALC_P_BONFERONI runs a Bonferoni test for significance in at least one
% period. The test is known to be very conservative.
% Source: Rom (1990)
% Authors: Makram El-Shagi, Axel Lindner and Gregor von Schweinitz
% _________________________________________________________________________
%
% SYNTAX:
%
% p = calc_p_Bonferoni(p_vec,alpha)
% _________________________________________________________________________
%
% INPUT
% 
% p_vec             Vector of observation-specific p-values
% alpha             significance level(s) to be tested
% _________________________________________________________________________
%
% OUTPUT
%
% p                 lowest alpha level where the test rejects

p_vec(p_vec==0) = eps;
p_vec(p_vec==1) = 1-eps;

N = size(p_vec,1);
if N==1
    p_vec = p_vec';
    N = size(p_vec,1);
end
p_vec2 = sort(p_vec,'descend');
reject = zeros(length(alpha),1);

logn = cumsum(log((1:N)));

for sigs = 1:length(alpha)
    alpha_N = zeros(N,1);
    alpha_N(1) = alpha(sigs);
    s1 = 0;
    for k = 2:N
        s1 = s1 + alpha(sigs)^(k-1);
        s2 = 0;
        if k>2
            for i = 1:k-2
%                 s2 = s2 + nchoosek(k,i)*alpha_N(i+1)^(k-i);
                add = (k-i)*log(alpha_N(i+1)) + logn(k) - logn(i) - logn(k-i);
                s2 = s2 + exp(add);
            end
        end
        alpha_N(k) = 1/k*(s1+s2);
    end
%     alpha_N(end-5:end)
    reject(sigs) = (sum(p_vec2<alpha_N)>0);
end

pos = find(reject,1);
if isempty(pos)
    p = 1;
else
    p = alpha(pos);
end