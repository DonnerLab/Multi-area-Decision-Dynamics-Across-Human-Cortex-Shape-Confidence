function T = recovery_tasklist(cfg)
%% =======================================================================
% One row per cluster task (one fit): [generating model, dataset, fitted
% model]. Results are stored per (g, d, f); a task whose result file
% exists is skipped.
%% =======================================================================

T = zeros(0,3);

for dataset = 1:cfg.n_datasets_model
    for fittedModel = cfg.models_in_recovery
        T(end+1,:) = [cfg.generating_model dataset fittedModel]; %#ok<AGROW>
    end
end

end
