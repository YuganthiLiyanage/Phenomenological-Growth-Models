clear all
close all
clc

   
numiter = 1000; 

     
% true_params = [0.000717,8.13,41.22,43,26.94]; 

true_params = [0.78, 0.85, 1.1e4]; 

        
 
X = zeros(length(true_params),numiter); 


dt = 0.1;

tforward = 0:dt:70;

t_measure =  (0:1:65)/dt + 1;
initial_cond = 3;


[~, y_trp] = ode15s(@(t,y)Model_GLM(y,true_params),tforward,initial_cond);

 Model_cases_C = y_trp(t_measure(:),1);
 Model_incidence =  0.78* (Model_cases_C.^0.85).*(1-Model_cases_C/1.1e4);


 noiselevel = [0.01, 0.05, 0.1, 0.2];

 total_ARE =  zeros(length(noiselevel), length(true_params));


 total_ARE_Table = {'r', 'alpha', 'k'};


for noisei = 1:4
    
rng default
noiselev = noiselevel(noisei)

    parfor i = 1:numiter
            i

  Incidence = (noiselev*randn(length(t_measure),1)) + Model_incidence;
  % Incidence = (noiselev*randn(length(t_measure),1)).*Model_incidence + Model_incidence;

   % LB=[rlb  0.01 1 20 I0 LBe];
   % UB=[rub  1 1 Kmax I0 UBe];

   % rlb=mean(abs(data1(1:2,2)))/200 = 0.0550;
   % rub=max(abs(data1(1:2,2)))*5=95;
 
            %params = [r alpha k]
            lb = [0.0550 0.01 20];  
            ub = [95 1 10000000000];
             
%              k = fmincon(@(k)err_in_data(k,Incidence), true_params,[], [], [], [], lb, ub,[],optimset('MaxFunEvals', 1e+5,'MaxIter',1e+5,'TolX',1e-12, 'Display','iter'))

             k = fmincon(@(k)err_in_data(k,Incidence), true_params,[], [], [], [], lb, ub,[],optimset('MaxFunEvals', 10000,'MaxIter',10000,'TolX',1.0000e-6));
             X(:,i) = k';
             
     end
        
        arescore = zeros(1,length(true_params));

    for i = 1:length(true_params)
        arescore(i) = 100*sum(abs(true_params(i) - X(i,:))/abs(true_params(i)))/numiter;
    end
    
    % total_ARE(noisei,:) = round(arescore,1);
    total_ARE(noisei,:) = arescore;
    total_ARE_Table(noisei+1,:) = num2cell(total_ARE(noisei,:));

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

tforward = 0:dt:70;

t_measure =  (0:1:65)/dt + 1;
initial_cond = 3;
 
 [~,y] = ode15s(@(t,y)Model_GLM(y,k),tforward,initial_cond);
 

 
 Model_C = y(t_measure(:),1);
 Model_Incidence_f=  k(1)* (Model_C.^k(2)).*(1-Model_C/k(3));

 
 error_in_data = sum((Model_Incidence_f - Incidence).^2);
                        

 end

 function dy = Model_GLM(y,k)

dy = zeros(1,1);

%params = [r alpha k]
r = k(1);
alpha = k(2);
d = k(3);




C = y(1);



dy(1) =  r* (C.^alpha).*(1-C/d);

 
end