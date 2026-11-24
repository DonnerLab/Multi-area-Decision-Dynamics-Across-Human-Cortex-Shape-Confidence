function [nLL, p, Prob_conf_choice, output] = cost_function_model(model_id, pars, data_in)
%COST_FUNCTION_MODEL  Common cost function for ALL models; calls model_<model_id>_simulations.m (unchanged).
%   Same likelihood as cost_function_model_9.m (probabilities smoothed as 0.99*P+0.005).
%
%   nLL = cost_function_model(9, pars)            % uses the global 'dat', exactly like cost_function_model_9
%   nLL = cost_function_model(9, pars, dat)       % explicit data struct (no globals)
%
%   dat fields: left, right (10 x nTrials), response (1 = right, 0 = left), conf (2 = high, 1 = low), t (optional)
%   Outputs (as in cost_function_model_9): p = P(choose right) per trial, Prob_conf_choice = n x 4
%   [R/high, R/low, L/high, L/low], output = P(high confidence) per trial.
%
%   Use in BADS:  bads(@(x) cost_function_model(9, x), x0, lb, ub, plb, pub, [], options)
global dat
if nargin < 3 || isempty(data_in)
    D = dat;
else
    D = data_in;
end
if ~isfield(D, 't'), D.t = size(D.left, 1); end

Prob_conf_choice = recovery_probs(model_id, pars, D, 1);
nLL    = -recovery_loglik(Prob_conf_choice, D.response, D.conf);
p      = Prob_conf_choice(:,1) + Prob_conf_choice(:,2);
output = (Prob_conf_choice(:,1) + Prob_conf_choice(:,3))';
end
