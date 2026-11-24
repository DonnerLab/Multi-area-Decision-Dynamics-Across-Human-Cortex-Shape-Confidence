function sig_mask = cluster_perm_1samp(X,nPerm,alpha)
% Cluster-based permutation test against zero (one sample)
%
% Inputs:
%   X     - time x subjects
%   nPerm - number of random sign flips of the subjects
%   alpha - cluster-forming threshold (|t| > tinv(1-alpha, nSubjects-1))
%           and cluster-level threshold ((1-alpha) percentile of the
%           maximum cluster mass under the null)
%
% Output:
%   sig_mask - time x 1, true = time point in a significant cluster

[nT,nS] = size(X);
t_obs = nan(nT,1);
for t=1:nT
    [~,~,~,stats] = ttest(X(t,:),0);
    t_obs(t) = stats.tstat;
end

thresh = tinv(1-alpha,nS-1);
clusters = bwlabel(abs(t_obs)>thresh);
clust_mass = accumarray(clusters(clusters>0),abs(t_obs(clusters>0)));

max_null = zeros(nPerm,1);
for p=1:nPerm
    flip = (rand(1,nS)>0.5)*2-1;
    t_null = nan(nT,1);
    for t=1:nT
        [~,~,~,stats] = ttest(X(t,:).*flip,0);
        t_null(t) = stats.tstat;
    end
    cl = bwlabel(abs(t_null)>thresh);
    if any(cl(:))
        max_null(p) = max(accumarray(cl(cl>0),abs(t_null(cl>0))));
    end
end

crit = prctile(max_null,100*(1-alpha));
sig_mask = false(nT,1);
for c=1:max(clusters)
    if clust_mass(c)>crit
        sig_mask(clusters==c)=true;
    end
end
end
