function ll = recovery_loglik(P, resp, conf)
%RECOVERY_LOGLIK  Log-likelihood of observed (choice, confidence) given model probabilities P (n x 4).
%   Identical to cost_function_model_9.m: probabilities smoothed as 0.99*P+0.005.
%   resp: 1 = right, 0 = left;  conf: 2 = high, 1 = low.
resp = resp(:);  conf = conf(:);
n    = numel(resp);
col  = 1 + 2*(resp==0) + (conf==1);        % 1: R/high  2: R/low  3: L/high  4: L/low
P    = max(P,0)*0.99 + 0.005;
lp   = log(P(sub2ind(size(P), (1:n)', col)));
ll   = sum(lp(~isnan(lp)));
end
