function model = simulate_model_dvs(parameters, data, simulationID, session)
%% =======================================================================
% Single-trial model DVs of one subject on the time base of the neural
% DVs (one simulation per trial)
%
% Inputs:
%   parameters   - model parameters (11 values)
%   data         - .left, .right: samples x trials contrasts (.t = 10)
%   simulationID - model variant of leaky_integrator_simulations_fine.m
%   session      - trials x 1, session of each trial
%
% Output (trials x 201 time points, 10 ms each, first at -0.5 s, as the
% decoder output; time point 51 = stimulus onset, 150 = 1 s):
%   model.DV_R, model.DV_L - right and left accumulator
%                            before stimulus onset: starting value
%                            (-bias for DV_R, +bias for DV_L)
%                            after stimulus offset: last value
%   model.isRight          - model choice (true = right)
%   model.conf             - model confidence (1 = high, 0 = low)
%   model.rt               - decision time (ms) = first bound crossing
%   model.ev               - evidence strength |mean(right - left)|
%   model.sess             - session of each trial
%
% Absorbing bound: the winning DV stays at the bound from its crossing
% on; the losing DV keeps accumulating until it reaches the bound itself
% (absorb_bound_dvs.m)
%% =======================================================================

nTrials     = size(data.left,2);
nTimePoints = 201;
onset       = 51;

[~, Choices, ~, modelConfidence, bound, ~, ~, P_left, P_right, ~, ~, Cross_step] = ...
    leaky_integrator_simulations_fine(parameters, data, simulationID, 1);

% Accumulators during the stimulus: time steps (10 ms) x trials
nSteps = size(P_left,1);

[P_left, P_right] = absorb_bound_dvs( ...
    reshape(P_left,nSteps,nTrials), ...
    reshape(P_right,nSteps,nTrials), ...
    bound, ...
    Choices(:) == 1, ...
    Cross_step(:));

stimulusTime = onset + (0:nSteps-1);
preStimulus  = 1:(onset-1);
postStimulus = (stimulusTime(end)+1):nTimePoints;

model.DV_R = zeros(nTrials,nTimePoints);
model.DV_L = zeros(nTrials,nTimePoints);

model.DV_R(:,preStimulus) = -parameters(11);
model.DV_L(:,preStimulus) =  parameters(11);

model.DV_R(:,stimulusTime) = P_right';
model.DV_L(:,stimulusTime) = P_left';

model.DV_R(:,postStimulus) = repmat(P_right(end,:)', 1, numel(postStimulus));
model.DV_L(:,postStimulus) = repmat(P_left(end,:)',  1, numel(postStimulus));

model.isRight = Choices(:) == 1;
model.conf    = double(modelConfidence(:));
model.rt      = (Cross_step(:) - 1) * 10;
model.ev      = abs(mean(data.right - data.left,1))';
model.sess    = session(:);

end
