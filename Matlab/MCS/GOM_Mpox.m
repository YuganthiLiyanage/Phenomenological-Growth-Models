clear all
close all
clc

   
numiter = 1000; 

     
% true_params = [0.000717,8.13,41.22,43,26.94]; 

true_params = [1.5 , 0.2]; 

        
 
X = zeros(length(true_params),numiter); 


dt = 0.1;

tforward = 0:dt:50;

t_measure =  (0:1:31)/dt + 1;
initial_cond = 14;


[t, y_trp] = ode15s(@(t,y)Model_GOM(t,y,true_params),tforward,initial_cond);

 Model_cases_C = (y_trp(t_measure(:),1))';
 Model_incidence = 1.5*Model_cases_C.*exp(-0.2*t_measure);

 noiselevel = [0.01, 0.05, 0.1, 0.2];

 total_ARE =  zeros(length(noiselevel), length(true_params));


 total_ARE_Table = {'r', 'b'};


for noisei = 1:4
    
rng default
noiselev = noiselevel(noisei)

    parfor i = 1:numiter
            i

  Incidence = noiselev*(randn(length(t_measure),1))' + Model_incidence;
  % Incidence = (noiselev*randn(length(t_measure),1)).*Model_incidence + Model_incidence;
  
  % LB=[0.0001  1 0 1 I0 LBe];
  % UB=[params0(1)+5  1 params0(3)+5 1 I0 UBe];
 
            
            lb = [0.0001 0];  
            ub = [6.5 5.2];
             
%              k = fmincon(@(k)err_in_data(k,Incidence), true_params,[], [], [], [], lb, ub,[],optimset('MaxFunEvals', 1e+5,'MaxIter',1e+5,'TolX',1e-12, 'Display','iter'))

             k = fmincon(@(k)err_in_data(k,Incidence), true_params,[], [], [], [], lb, ub,[],optimset('MaxFunEvals', 10000,'MaxIter',10000,'TolX',1.0000e-6));
             X(:,i) = k';
             
     end
        
        arescore = zeros(1,length(true_params));

    for i = 1:length(true_params)
        arescore(i) = 100*sum(abs(true_params(i) - X(i,:))/abs(true_params(i)))/numiter;
    end
    
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

tforward = 0:dt:50;

t_measure =  (0:1:31)/dt + 1;
initial_cond = 14;
 
 [~,y] = ode15s(@(t,y)Model_GOM(t,y,k),tforward,initial_cond);
 

 
 Model_C = y(t_measure(:),1);
 Model_Incidence_f= k(1)* (Model_C)'.*exp(-k(2)*t_measure);
 
 error_in_data = sum((Model_Incidence_f - Incidence).^2);
                        

 end

 function dy = Model_GOM(t,y,k)

 dy = zeros(1,1);


%params = [r b]
r = k(1);
b = k(2);

% dt=0.1;
% tforward = 0:dt:50;
% t_measure =  (0:1:31)/dt + 1;


C = y(1);



dy(1) =  r*C* exp(-b*t);

 
end