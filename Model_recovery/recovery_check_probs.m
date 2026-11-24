function recovery_check_probs(P, m)
%RECOVERY_CHECK_PROBS  Stop early if model m returns invalid category probabilities.
if size(P,2) ~= 4 || any(isnan(P(:))) || any(abs(sum(P,2) - 1) > 1e-3)
    error(['model_%d_simulations returned invalid probabilities (rows must sum to 1). ' ...
           'For the no-bound models (1,4,7,12,15) check that every simulation gets a choice ' ...
           '(Ch must be nan(size(P)), not ones(size(P))+9).'], m);
end
end
