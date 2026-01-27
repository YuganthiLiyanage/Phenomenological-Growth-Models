clear all
close all
clc

   
numiter = 1000; 


% parms=['beta','k','gamma','N', 'alpha']; 

true_params = [7.3,4.7,4.8,1e5,0.96]; 
% true_params = [6.5,0.5,0.33,1e5,0.001]; 

        
X = zeros(length(true_params),numiter); 


dt = 0.1;

tforward = 0:dt:32;

t_measure =  (0:1:31)/dt + 1;
initial_cond = [99995 0 5 0];


[~,y_trp] = ode15s(@(t,y)Model_SEIR(y,true_params), tforward,initial_cond);

 Model_E = y_trp(t_measure(:),2);
 Model_Obs = true_params(2)*Model_E;

 noiselevel = [0, 0.01, 0.05, 0.1, 0.2];

 total_ARE =  zeros(length(noiselevel), length(true_params));


 total_ARE_Table = {'beta','k','gamma','N', 'alpha'};


for noisei = 1:5
    
rng default
noiselev = noiselevel(noisei)

    parfor i = 1:numiter
            i

  Sim_Data = noiselev*randn(length(t_measure),1) + Model_Obs;
  
 
          % parms=['beta','k','gamma','N', 'alpha'];  
             lb = [0.01 0.01 0.01 20 0];  
            ub = [25 5 5 1000000 1]; 
             
           
              k = fmincon(@(k)err_in_data(k,Sim_Data), true_params,[], [], [], [], lb, ub,[],optimset('MaxFunEvals', 10000,'MaxIter',10000,'TolX',1.0000e-6));
             X(:,i) = k';
             
     end
        
        arescore = zeros(1,length(true_params));

    for i = 1:length(true_params)
        arescore(i) = 100*sum(abs(true_params(i) - X(i,:))/abs(true_params(i)))/numiter;
    end
    
    total_ARE(noisei,:) = arescore;
    total_ARE_Table(noisei+1,:) = num2cell(total_ARE(noisei,:));

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

function error_in_data = err_in_data(k,Sim_Data) 
 
dt = 0.1;

tforward = 0:dt:32;

t_measure =  (0:1:31)/dt + 1;
initial_cond = [99995 0 5 0];



 
 [~,y] = ode15s(@(t,y)Model_SEIR(y,k), tforward,initial_cond);
 

 
 Model_Obs_f =k(2)* y(t_measure(:),2);


  error_in_data=sum((Model_Obs_f - Sim_Data).^2);
 
 

                        

 end

 function dy = Model_SEIR(y,k)

dy = zeros(4,1);

% parms=['beta','k','gamma','N', 'alpha']; 
beta = k(1);
d = k(2);
gamma = k(3);
N = k(4);
alpha = k(5);



S = y(1);
E = y(2);
I = y(3);
R = y(4);


dy(1) = - (beta* S.*I.^alpha)./N ;
dy(2) = (beta* S.*I.^alpha)./N  - d*E;
dy(3) = d*E - gamma*I;
dy(4) = gamma*I;
 
end