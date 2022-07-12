function err_norm2 = cyclic_error_calculation(optimized_mat_paras)

% This function is used to calculate the total differnence betweeen the test
% results and model prediction of all required coupons

%% set global variables
global optimization_data_name_set ... % coupon name list for optimization
       elastic_paras ... % elastic parameters, i.e., elastic modulus
       monotonic_paras ... % material parameters for monotonic loading
       elastic_modulus_paras ... % parameters for elastic modulus evolution

%% generate material prameter vaector   
mat_paras = [elastic_paras, monotonic_paras, optimized_mat_paras, elastic_modulus_paras];

%% calculate error of each test
% intialize the error list
num_optim_dataset = numel(optimization_data_name_set);
error_list = zeros(num_optim_dataset,1);

% compute error
for i  = 1:num_optim_dataset

    name_of_loading = optimization_data_name_set{i};
    error_list(i) = cyclic_errori(mat_paras,name_of_loading);

end

%% calculate the total error
err_norm2 = sum(error_list);

end
