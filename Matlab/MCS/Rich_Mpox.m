clear all
close all
clc

numiter = 1000;  

% [r, a, k]
true_params = [0.92, 0.36, 2.9e4];  


X = zeros(length(true_params), numiter);


dt = 0.1;
tforward = 0:dt:50;
t_measure = (0:1:31)/dt + 1;
initial_cond = 14;


[~, y_trp] = ode15s(@(t,y)Model_Rich(y,true_params), tforward, initial_cond);
Model_cases_C = y_trp(t_measure(:),1);
Model_incidence = 0.92 * Model_cases_C .* (1 - (Model_cases_C/2.9e4).^0.36);


noiselevel = [0, 0.01, 0.05, 0.1, 0.2];
total_ARE = zeros(length(noiselevel), length(true_params));
total_ARE_Table = {'r', 'a', 'k'};


all_param_medians = zeros(length(true_params), length(noiselevel));
all_conf_intervals_lower = zeros(length(true_params), length(noiselevel));
all_conf_intervals_upper = zeros(length(true_params), length(noiselevel));


for noisei = 1:length(noiselevel)
    rng default
    noiselev = noiselevel(noisei);

    parfor i = 1:numiter
        i
        
        % Add noise
        Incidence = (noiselev * randn(length(t_measure), 1)) + Model_incidence;
        
        
        lb = [0.0750 0 0];  
        ub = [80 10 1e10];
        
       
        k = fmincon(@(k)err_in_data(k, Incidence), true_params, [], [], [], [], lb, ub, [], ...
                    optimset('MaxFunEvals', 10000, 'MaxIter', 10000, 'TolX', 1.0000e-6, 'Display', 'off'));
        
        
        X(:, i) = k';
    end
    
   
    arescore = zeros(1, length(true_params));
    for j = 1:length(true_params)
        arescore(j) = 100 * sum(abs(true_params(j) - X(j,:)) / abs(true_params(j))) / numiter;
    end
    total_ARE(noisei,:) = arescore;
    total_ARE_Table(noisei + 1,:) = num2cell(total_ARE(noisei,:));

    % Calculate confidence intervals (95%) 
    param_medians = median(X, 2);  % Median 
    conf_intervals_lower = quantile(X, 0.025, 2);  % 2.5th percentile 
    conf_intervals_upper = quantile(X, 0.975, 2);  % 97.5th percentile 

    % Store 
    all_param_medians(:, noisei) = param_medians;
    all_conf_intervals_lower(:, noisei) = conf_intervals_lower;
    all_conf_intervals_upper(:, noisei) = conf_intervals_upper;
end


fprintf('Final Results:\n');
for noisei = 1:length(noiselevel)
    fprintf('\nNoise Level = %.2f:\n', noiselevel(noisei));
    for j = 1:length(true_params)
        fprintf('Parameter %d: Median = %.4f, 95%% Confidence Interval: [%.4f, %.4f]\n', ...
                j, all_param_medians(j, noisei), all_conf_intervals_lower(j, noisei), all_conf_intervals_upper(j, noisei));
    end
end

% Functions
function error_in_data = err_in_data(k, Incidence)
    dt = 0.1;
    tforward = 0:dt:50;
    t_measure = (0:1:31)/dt + 1;
    initial_cond = 14;
    
    [~, y] = ode15s(@(t,y)Model_Rich(y,k), tforward, initial_cond);
    Model_C = y(t_measure(:),1);
    Model_Incidence_f = k(1) * Model_C .* (1 - (Model_C/k(3)).^k(2));
    
    error_in_data = sum((Model_Incidence_f - Incidence).^2);
end

function dy = Model_Rich(y, k)
    dy = zeros(1,1);
    
    r = k(1);
    a = k(2);
    d = k(3);
    
    C = y(1);
    dy(1) = r * C .* (1 - (C/d).^a);
end
