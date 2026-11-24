function [Final_Confidence, Choices, Confidence, output, bound, ...
          Confidence_left, Confidence_right, P_left, P_right, ...
          P_cross_id, Cross_time, Cross_step] = ...
          leaky_integrator_simulations_fine(pars, data, m_n, n_simulations)
%% =======================================================================
% Accumulator model simulated with a time step 10 times finer than the
% stimulus samples (100 steps of 10 ms during the 1 s stimulus), with the
% parameters fitted with the 10-step model (model_fit.m)
%
% Refinement of the 10-step model:
%   1. adaptation evaluated continuously (at each fine step)
%   2. collapsing bound evaluated continuously
%   3. each contrast sample held constant for 10 fine steps
%   4. leak rescaled so that 10 fine steps retain as much as one sample
%   5. noise rescaled so that the noise variance accumulated over 10 fine
%      steps equals the noise variance of one sample
%   6. bound crossing detected at each fine step
%
% Inputs:
%   pars          - fitted model parameters
%   data          - .left, .right: samples x trials contrasts
%   m_n           - model variant:
%                     2 = fixed bound, feed-forward inhibition (FFI)
%                     4 = collapsing bound, FFI, confidence asymmetry
%                         weight (best-fitting model, model 9)
%                     5 = collapsing bound, FFI, rectification of the
%                         accumulators at 0 (model 14)
%   n_simulations - number of simulations per trial
%
% Outputs:
%   Choices          - 1 x trials x simulations (0 = left, 1 = right)
%   output           - 1 x trials x simulations (0 = low, 1 = high
%                      confidence)
%   Confidence       - trials x 4 probabilities:
%                      [right/high, right/low, left/high, left/low]
%   bound            - bound at each fine step
%   P_left, P_right  - fine steps x trials x simulations, accumulators
%   Cross_time       - 1 x trials x simulations, time of the first bound
%                      crossing in units of stimulus samples (0 ... 9.9);
%                      NaN = no crossing
%   Cross_step       - fine step of the first bound crossing (NaN = no
%                      crossing)
%   P_cross_id       - fine steps from the crossing on (1), otherwise 0
%   Final_Confidence - |confidence variable| before the logistic readout
%   Confidence_left, Confidence_right - accumulators at the readout
%                      (winning accumulator = bound if it crossed)
%% =======================================================================


%% =======================================================================
% PARAMETERS
%% =======================================================================

dt      = 0.1;    % fine step, in units of stimulus samples
nSteps  = 100;
time    = (0:nSteps-1) * dt;

tau_a   = pars(3) * 15 + 0.0001;   % adaptation time constant
noise   = pars(2) / 20;            % noise SD per stimulus sample
leak    = pars(1);
w_I     = pars(5);                 % feed-forward inhibition
bias    = pars(11);

% Confidence asymmetry weight (variant 4 only)
if m_n == 4
    w_conf = pars(10);
else
    w_conf = 1;
end


%% =======================================================================
% LEAK AND NOISE OF ONE FINE STEP
%
% Retention over one sample: 1 - leak = a^10, with a the retention over
% one fine step. Noise SD of one fine step chosen so that the variance
% accumulated over 10 fine steps equals noise^2:
%   noise_fine^2 * sum(a^(2j), j = 0:9) = noise^2
%% =======================================================================

a_fine = exp(log(1 - leak) * dt);

if abs(1 - a_fine^2) < eps
    noise_fine = noise * sqrt(1/10);
else
    noise_fine = noise * sqrt((1 - a_fine^2) / (1 - a_fine^20));
end


%% =======================================================================
% ADAPTATION AND BOUND AT EACH FINE STEP
%% =======================================================================

% Adaptation kernel (normalized to sum 1 over the stimulus)
adaptation = pars(4) + (1 - pars(4)) .* exp(-time ./ tau_a);
adaptation = adaptation ./ sum(adaptation);

% Bound: fixed (variant 2) or collapsing
initialBound = abs(pars(8)) * 15;
collapseRate = abs(pars(6) * 50);

if m_n == 2
    bound = ones(1,nSteps) * initialBound;
else
    bound = initialBound - initialBound .* (time ./ (time + 1/collapseRate));
end


%% =======================================================================
% INPUT: each contrast sample held for 10 fine steps, times adaptation
%% =======================================================================

nTrials = size(data.left,2);

left_input  = repelem(data.left(1:10,:),  10, 1) .* adaptation';
right_input = repelem(data.right(1:10,:), 10, 1) .* adaptation';


%% =======================================================================
% ACCUMULATION
%
% P(k) = a_fine * P(k-1) + input(k) - w_I * input of the other side(k)
%        + noise; variant 5: accumulators rectified at 0
%% =======================================================================

P_left  = zeros(nSteps, nTrials, n_simulations);
P_right = zeros(nSteps, nTrials, n_simulations);

for k = 1:nSteps

    noise_left  = normrnd(0, noise_fine, [1 nTrials n_simulations]);
    noise_right = normrnd(0, noise_fine, [1 nTrials n_simulations]);

    if k == 1
        previous_left  = 0;
        previous_right = 0;
    else
        previous_left  = a_fine .* P_left(k-1,:,:);
        previous_right = a_fine .* P_right(k-1,:,:);
    end

    P_left(k,:,:)  = previous_left  + (left_input(k,:)  - w_I .* right_input(k,:)) + noise_left;
    P_right(k,:,:) = previous_right + (right_input(k,:) - w_I .* left_input(k,:))  + noise_right;

    if m_n == 5
        P_left(P_left < 0)   = 0;
        P_right(P_right < 0) = 0;
    end

end

% Bias
P_left  = P_left  + bias;
P_right = P_right - bias;


%% =======================================================================
% FIRST BOUND CROSSING
%
% Readout at the first crossing of either accumulator, or at the last
% fine step if the bound is never reached
%% =======================================================================

crossed_left  = P_left  >= bound';
crossed_right = P_right >= bound';
crossed_any   = crossed_left | crossed_right;

has_crossing = any(crossed_any,1);
[~, first_step] = max(crossed_any,[],1);
first_step(~has_crossing) = nSteps;

% Values at the readout step of every trial x simulation
readout = first_step + reshape((0:(nTrials*n_simulations-1)) .* nSteps, 1, nTrials, n_simulations);

P_left_read   = P_left(readout);
P_right_read  = P_right(readout);
crossed_left  = crossed_left(readout);
crossed_right = crossed_right(readout);
bound_read    = reshape(bound(first_step), 1, nTrials, n_simulations);


%% =======================================================================
% CHOICE (+1 = right, -1 = left)
%
% One accumulator crossed: its side; both crossed: the larger one (ties:
% left); no crossing: sign of right - left (ties: random)
%% =======================================================================

choice = nan(1,nTrials,n_simulations);

choice(has_crossing & crossed_right & ~crossed_left) =  1;
choice(has_crossing & crossed_left  & ~crossed_right) = -1;

both = has_crossing & crossed_left & crossed_right;
choice(both & P_right_read >  P_left_read) =  1;
choice(both & P_right_read <= P_left_read) = -1;

choice(~has_crossing) = sign(P_right_read(~has_crossing) - P_left_read(~has_crossing));

ties = find(choice == 0);
choice(ties) = sign(randn(size(ties)));


%% =======================================================================
% CROSSING TIME
%% =======================================================================

Cross_step = nan(1,nTrials,n_simulations);
Cross_time = nan(1,nTrials,n_simulations);

Cross_step(has_crossing) = first_step(has_crossing);
Cross_time(has_crossing) = time(first_step(has_crossing));

P_cross_id = double((reshape(1:nSteps,[],1,1) >= first_step) & has_crossing);


%% =======================================================================
% CONFIDENCE VARIABLE
%
% Crossing   : bound - w_conf * losing accumulator
% No crossing: |winning - w_conf * losing accumulator| at the last step
%% =======================================================================

confidenceVariable = nan(1,nTrials,n_simulations);
Confidence_left    = P_left_read;
Confidence_right   = P_right_read;

right = choice ==  1;
left  = choice == -1;

id = ~has_crossing & right;
confidenceVariable(id) = abs(P_right_read(id) - w_conf .* P_left_read(id));

id = ~has_crossing & left;
confidenceVariable(id) = abs(P_left_read(id) - w_conf .* P_right_read(id));

id = has_crossing & right;
confidenceVariable(id) = bound_read(id) - w_conf .* P_left_read(id);
Confidence_right(id)   = bound_read(id);

id = has_crossing & left;
confidenceVariable(id) = bound_read(id) - w_conf .* P_right_read(id);
Confidence_left(id)    = bound_read(id);

% Choice: 0 = left, 1 = right
choice(choice == -1) = 0;


%% =======================================================================
% CONFIDENCE REPORT
%
% |confidence variable| normalized between its minimum and maximum
% across trials (separately for each simulation), then logistic function
% with slope pars(7) and threshold pars(9) x median across trials:
% high confidence if > 0.5 (exactly 0.5: random)
%% =======================================================================

slope_conf = abs(pars(7) * 30);
thr_conf   = abs(pars(9) * 4);

Final_Confidence = abs(confidenceVariable);

minimum = min(Final_Confidence,[],2);
maximum = max(Final_Confidence,[],2);

normalizedConfidence = (Final_Confidence - minimum) ./ (maximum - minimum);

pHigh = 1 ./ (1 + exp(-slope_conf .* (normalizedConfidence - nanmedian(normalizedConfidence,2) .* thr_conf)));

exactlyHalf = find(abs(pHigh) == 0.5);
output = double(pHigh > 0.5);
output(exactlyHalf) = randi(2, size(exactlyHalf)) - 1;


%% =======================================================================
% PROBABILITIES OF THE 4 REPORTS (across simulations)
%% =======================================================================

Confidence(:,1) = sum(choice == 1 & output == 1, 3) ./ n_simulations;
Confidence(:,2) = sum(choice == 1 & output == 0, 3) ./ n_simulations;
Confidence(:,3) = sum(choice == 0 & output == 1, 3) ./ n_simulations;
Confidence(:,4) = sum(choice == 0 & output == 0, 3) ./ n_simulations;

Choices = choice;

end
