function [nLL,p,Prob_conf_choice,output] = cost_function_model_9(pars)

%% Cost function model 9

global dat

% Real responses
resp      = dat.response;
data_conf = dat.conf;

% Calculate model predictions
[p, Prob_conf_choice, output] = model_9_simulations(pars,dat);

% Matrix of observed choice-confidence combinations:
%   col 1 = choice 1, high confidence
%   col 2 = choice 1, low confidence
%   col 3 = choice 0, high confidence
%   col 4 = choice 0, low confidence
Data_choice = zeros(size(p,1),4);

for tr = 1:size(p,1)

    if resp(tr)==1 && data_conf(tr)==2
        Data_choice(tr,1)=1;

    elseif resp(tr)==1 && data_conf(tr)==1
        Data_choice(tr,2)=1;

    elseif resp(tr)==0 && data_conf(tr)==2
        Data_choice(tr,3)=1;

    else
        Data_choice(tr,4)=1;
    end

end

% Avoid zero probabilities
Prob_conf_choice(Prob_conf_choice<0) = 0;

Choices_P = double(Prob_conf_choice)*0.99 + 0.005;

% Log likelihood
prob_choice = log(nansum(Data_choice .* Choices_P,2));

nLL = -nansum(prob_choice);

end