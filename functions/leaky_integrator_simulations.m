
function [Final_Confidence, Choices, Confidence, output, bound, ...
          Confidence_left, Confidence_right, P_left, P_right, P_cross_id] = ...
          leaky_integrator_simulations(pars, data, m_n, n_simulations)

%% ========================================================================
% Leaky integrator model simulations
%
% Model variants (m_n):
% 1 = No bound, FFI, no confidence asymmetry
% 2 = Fixed bound, FFI, no confidence asymmetry
% 3 = Collapsing bound, FFI, no confidence asymmetry
% 4 = Collapsing bound, FFI, confidence asymmetry
% 5 = Collapsing bound, FFI, rectification at 0
% 6 = No bound, LCA, rectification at 0
% 7 = Collapsing bound, FFI, reflective bound
%% ========================================================================

%% ========================================================================
% Adaptation kernel
%% ========================================================================

tau_a = pars(3) * 15 + 0.0001;

Adapt = pars(4) + (1 - pars(4)) * exp(-(0:9) ./ tau_a);
Adapt = Adapt ./ sum(Adapt);

%% ========================================================================
% Basic parameters
%% ========================================================================

t    = data.t - 1;
ext  = pars(2) / 20;      % external noise
leak = pars(1);           % leak parameter
w_I  = pars(5);           % inhibition weight
B    = pars(11);          % starting point bias

%% ========================================================================
% Confidence asymmetry weight
%% ========================================================================

if m_n == 4
    w_conf = pars(10);
else
    w_conf = 1;
end

%% ========================================================================
% Decision bound
%% ========================================================================

dv    = abs(pars(6) * 50);
thr_1 = abs(pars(8)) * 15;

if ismember(m_n, [1 6])
    
    % No bound
    bound = inf(1,100) * thr_1;

elseif m_n == 2
    
    % Fixed bound
    bound = ones(1,100) * thr_1;

else
    
    % Collapsing bound
    time_axis = 0:0.1:9.9;
    bound = thr_1 - (thr_1 - 0) .* ...
        (time_axis ./ (time_axis + (1 ./ dv)));
end

%% ========================================================================
% Stimulus preprocessing + interpolation
%% ========================================================================

time_points     = 1:10;
new_time_points = linspace(1,10,100);

left_input = interp1( ...
    time_points, ...
    data.left(1:10,:) .* Adapt', ...
    new_time_points, ...
    'linear') * 0.1;

right_input = interp1( ...
    time_points, ...
    data.right(1:10,:) .* Adapt', ...
    new_time_points, ...
    'linear') * 0.1;

%% ========================================================================
% Initialize decision variables
%% ========================================================================

n_trials = size(left_input,2);

P_left  = zeros(100, n_trials, n_simulations);
P_right = zeros(100, n_trials, n_simulations);

%% ========================================================================
% First sample
%% ========================================================================

k = 1;

Noise_L = normrnd(0, ext/sqrt(10), [1 n_trials n_simulations]);
Noise_R = normrnd(0, ext/sqrt(10), [1 n_trials n_simulations]);

if m_n == 6
    % LCA
    L_to_R = 0;
    R_to_L = 0;
else
    % Feedforward inhibition
    L_to_R = w_I .* left_input(k,:,:);
    R_to_L = w_I .* right_input(k,:,:);
end

P_left(k,:,:)  = (left_input(k,:,:)  - R_to_L) + Noise_L;
P_right(k,:,:) = (right_input(k,:,:) - L_to_R) + Noise_R;

%% ========================================================================
% Rectification rules
%% ========================================================================

if ismember(m_n, [5 6])
    
    % Zero rectification
    P_left(P_left < 0)   = 0;
    P_right(P_right < 0) = 0;

elseif m_n == 7
    
    % Reflective lower bound
    P_left(P_left < pars(12))   = pars(12);
    P_right(P_right < pars(12)) = pars(12);
end

%% ========================================================================
% Evidence accumulation
%% ========================================================================

for k = 2:100

    Noise_L = normrnd(0, ext/sqrt(10), [1 n_trials n_simulations]);
    Noise_R = normrnd(0, ext/sqrt(10), [1 n_trials n_simulations]);

    if m_n == 6
        
        % LCA: inhibition depends on previous DV
        L_to_R = w_I .* P_left(k-1,:,:);
        R_to_L = w_I .* P_right(k-1,:,:);

    else
        
        % FFI: inhibition depends on current sensory input
        L_to_R = w_I .* left_input(k,:,:);
        R_to_L = w_I .* right_input(k,:,:);
    end

    P_left(k,:,:) = ...
        (1 - leak/10) .* P_left(k-1,:,:) + ...
        (left_input(k,:,:) - R_to_L) + Noise_L;

    P_right(k,:,:) = ...
        (1 - leak/10) .* P_right(k-1,:,:) + ...
        (right_input(k,:,:) - L_to_R) + Noise_R;

    % Rectification
    if ismember(m_n, [5 6])

        P_left(P_left < 0)   = 0;
        P_right(P_right < 0) = 0;

    elseif m_n == 7

        P_left(P_left < pars(12))   = pars(12);
        P_right(P_right < pars(12)) = pars(12);
    end
end

%% ========================================================================
% Bias
%% ========================================================================

P_left  = P_left  + B;
P_right = P_right - B;

%% ========================================================================
% Bound crossing
%% ========================================================================

Ch_end = nan(1, n_trials, n_simulations);

Ch_L = nan(size(P_left));
Ch_R = nan(size(P_right));

Ch_L(P_left >= bound')  = -1;
Ch_R(P_right >= bound') = 1;

id_both_nan = isnan(Ch_L) & isnan(Ch_R);

Ch = nansum(cat(4, Ch_L, Ch_R), 4);
Ch(id_both_nan) = nan;

%% ========================================================================
% Confidence containers
%% ========================================================================

Confidence_end   = nan(1, n_trials, n_simulations);
Confidence_left  = nan(1, n_trials, n_simulations);
Confidence_right = nan(1, n_trials, n_simulations);

P_cross    = zeros(1, n_trials, n_simulations);
P_cross_id = zeros(size(Ch));

%% ========================================================================
% Trials without bound crossing
%% ========================================================================

for tr = 1:n_trials

    final_sample = 100;

    id_end = find(sum(isnan(Ch(:,tr,:)),1) == 100);

    Ch_end(1,tr,id_end) = sign( ...
        P_right(final_sample,tr,id_end) - ...
        P_left(final_sample,tr,id_end));

    P_cross(1,tr,id_end) = 1;

    % Random tie-break
    id_zero = find(Ch_end(1,tr,id_end) == 0);
    Ch_end(1,tr,id_end(id_zero)) = sign(randn(size(id_zero)));

    % Confidence
    id_pos = find(Ch_end(1,tr,id_end) == 1);
    id_neg = find(Ch_end(1,tr,id_end) == -1);

    Confidence_end(1,tr,id_end(id_pos)) = ...
        abs(P_right(final_sample,tr,id_end(id_pos)) - ...
        w_conf .* P_left(final_sample,tr,id_end(id_pos)));

    Confidence_end(1,tr,id_end(id_neg)) = ...
        abs(P_left(final_sample,tr,id_end(id_neg)) - ...
        w_conf .* P_right(final_sample,tr,id_end(id_neg)));

    Confidence_left(1,tr,id_end)  = P_left(final_sample,tr,id_end);
    Confidence_right(1,tr,id_end) = P_right(final_sample,tr,id_end);
end

P_cross = 1 - P_cross;

%% ========================================================================
% Trials with bound crossing
%% ========================================================================

for i = 1:100

    this_Choice  = Ch(i,:,:);
    this_P_right = P_right(i,:,:);
    this_P_left  = P_left(i,:,:);

    id = find(this_Choice == 1 | this_Choice == -1);

    Ch_end(id) = this_Choice(id);

    id_equal = find(this_Choice == 0);

    Ch_end(id_equal) = double( ...
        this_P_right(id_equal) > this_P_left(id_equal));

    this_Choice(id_equal) = Ch_end(id_equal);
    this_Choice(this_Choice == 0) = -1;

    if m_n > 1
        this_P_right(id_equal) = bound(i);
        this_P_left(id_equal)  = bound(i);
    end

    Ch_end(Ch_end == 0) = -1;

    % Confidence
    id_pos = find(this_Choice == 1);
    id_neg = find(this_Choice == -1);

    Confidence_end(id_pos) = ...
        bound(i) - w_conf .* this_P_left(id_pos);

    Confidence_end(id_neg) = ...
        bound(i) - w_conf .* this_P_right(id_neg);

    Confidence_left(id_pos)  = this_P_left(id_pos);
    Confidence_right(id_neg) = this_P_right(id_neg);

    Confidence_left(id_neg)  = bound(i);
    Confidence_right(id_pos) = bound(i);

    % Mark bound crossing
    id_cross = (Ch(i,:,:) == 1 | Ch(i,:,:) == -1 | Ch(i,:,:) == 0);

    P_cross_id(i:end,id_cross) = 1;
    Ch(i+1:end,id_cross) = nan;
end

%% ========================================================================
% Final confidence report
%% ========================================================================

Ch_end(Ch_end == -1) = 0;

slope_conf = abs(pars(7) * 30);
thr_conf   = abs(pars(9) * 4);

Conf_all = abs(Confidence_end);
Final_Confidence = Conf_all;

min_v = min(min(Conf_all,[],2),[],1);
max_v = max(max(Conf_all,[],2),[],1);

min_v = repmat(min_v,1,size(Conf_all,2),1);
max_v = repmat(max_v,1,size(Conf_all,2),1);

Norm_Conf_all = (Conf_all - min_v) ./ (max_v - min_v);

Median = nanmedian(Norm_Conf_all,2);

Confidence_end = 1 ./ (1 + exp( ...
    -slope_conf .* (Norm_Conf_all - Median .* thr_conf)));

id_equal = find(abs(Confidence_end) == 0.5);

Confidence_end = double(Confidence_end > 0.5);
Confidence_end(id_equal) = randi(2,size(id_equal)) - 1;

%% ========================================================================
% Confidence probabilities
%% ========================================================================

id_1_1 = (Ch_end(1,:,:) == 1 & Confidence_end(1,:,:) == 1);
id_1_0 = (Ch_end(1,:,:) == 1 & Confidence_end(1,:,:) == 0);
id_0_1 = (Ch_end(1,:,:) == 0 & Confidence_end(1,:,:) == 1);
id_0_0 = (Ch_end(1,:,:) == 0 & Confidence_end(1,:,:) == 0);

P_confidence(:,1) = sum(id_1_1,3) ./ n_simulations;
P_confidence(:,2) = sum(id_1_0,3) ./ n_simulations;
P_confidence(:,3) = sum(id_0_1,3) ./ n_simulations;
P_confidence(:,4) = sum(id_0_0,3) ./ n_simulations;

%% ========================================================================
% Outputs
%% ========================================================================

Choices    = Ch_end;
output     = Confidence_end;
Confidence = P_confidence;

end