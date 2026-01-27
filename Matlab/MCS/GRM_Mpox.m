clear all
close all
clc

   
numiter = 1000; 

     
% true_params = [0.000717,8.13,41.22,43,26.94]; 

true_params = [2.3, 0.82, 3.1e4, 0.9]; 

        
 
X = zeros(length(true_params),numiter); 


dt = 0.1;

tforward = 0:dt:50;

t_measure =  (0:1:31)/dt + 1;
initial_cond = 14;


[~, y_trp] = ode15s(@(t,y)Model_GRM(y,true_params),tforward,initial_cond);

 Model_cases_C = y_trp(t_measure(:),1);
 Model_incidence = 2.3* (Model_cases_C.^0.82).*(1-(Model_cases_C/3.1e4).^0.9);

 noiselevel = [0,0.01, 0.05, 0.1, 0.2];

 total_ARE =  zeros(length(noiselevel), length(true_params));


 total_ARE_Table = {'r', 'alpha', 'k', 'a'};


for noisei = 1:5
    
rng default
noiselev = noiselevel(noisei)

    parfor i = 1:numiter
            i

  Incidence = (noiselev*randn(length(t_measure),1)) + Model_incidence;
  % Incidence = (noiselev*randn(length(t_measure),1)).*Model_incidence + Model_incidence;
  
 % LB=[rlb  0.01 0 20 I0 LBe];
 % UB=[rub  1 10 Kmax I0 UBe];
% rlb=mean(abs(data1(1:2,2)))/200 = 0.0750;
% rub=max(abs(data1(1:2,2)))*5 = 80;
% Kmax=10000000000;
            %params = [r alpha k a]
            lb = [0.0750 0.01 20 0];  
            ub = [80 1 10000000000 10];
             
              k = fmincon(@(k)err_in_data(k,Incidence), true_params,[], [], [], [], lb, ub,[],optimset('MaxFunEvals', 10000,'MaxIter',10000,'TolX',1.0000e-6))
             % options=optimoptions('fmincon','Algorithm','sqp','StepTolerance',1.0000e-6,'MaxFunEvals',20000,'MaxIter',20000);
             % k = fmincon(@(k)err_in_data(k,Incidence), true_params,[], [], [], [], lb, ub,options)
             X(:,i) = k';
             
     end
        
        arescore = zeros(1,length(true_params));

    for i = 1:length(true_params)
        arescore(i) = 100*sum(abs(true_params(i) - X(i,:))/abs(true_params(i)))/numiter;
    end
    
    % total_ARE(noisei,:) = round(arescore,1);
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

function error_in_data = err_in_data(k, Incidence) 
 
dt = 0.1;

tforward = 0:dt:50;

t_measure =  (0:1:31)/dt + 1;
initial_cond = 14;
 
 [~,y] = ode15s(@(t,y)Model_GRM(y,k),tforward,initial_cond);
 

 
 Model_C = y(t_measure(:),1);
 Model_Incidence_f= k(1)* (Model_C.^k(2)).*(1-(Model_C/k(3)).^k(4));
 
 error_in_data = sum((Model_Incidence_f - Incidence).^2);
                        

 end

 function dy = Model_GRM(y,k)

dy = zeros(1,1);

%params = [r alphs k a]
r = k(1);
alpha = k(2);
d = k(3);
a = k(4);




C = y(1);



dy(1) =  r* (C.^alpha).*(1-(C/d).^a);

 
end