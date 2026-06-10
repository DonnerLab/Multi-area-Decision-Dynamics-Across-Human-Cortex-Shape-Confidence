function [h, qvals] = fdr_bh(p, q)
    % Benjamini-Hochberg FDR
    %
    % Inputs:
    %   p = vector of p-values
    %   q = FDR level (e.g. 0.01)
    %
    % Outputs:
    %   h = 1 if significant, 0 otherwise
    %   qvals = FDR adjusted p-values

    p = p(:);
    V = numel(p);

    [ps, idx] = sort(p);
    ranks = (1:V)';

    % BH critical values
    crit = (ranks/V) * q;

    % Find largest passing index
    below = find(ps <= crit);
    if isempty(below)
        h_sorted = false(V,1);
    else
        h_sorted = ranks <= max(below);
    end

    % Reorder h to original
    h = false(V,1);
    h(idx) = h_sorted;

    % Adjusted p-values (q-values)
    qvals_sorted = ps .* V ./ ranks;
    qvals_sorted(qvals_sorted > 1) = 1;

    % Enforce monotonicity (required)
    qvals_sorted = cummin(qvals_sorted(end:-1:1));
    qvals_sorted = qvals_sorted(end:-1:1);

    % Return qvals in original order
    qvals = zeros(V,1);
    qvals(idx) = qvals_sorted;

    % Reshape
    h = reshape(h, size(p));
    qvals = reshape(qvals, size(p));
end
