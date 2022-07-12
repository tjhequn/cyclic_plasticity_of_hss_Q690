function err_norm2 = monotonic_error_calculation(optimized_mono_mat_paras)

% This function is used to calculate the total differnence betweeen the test
% results and model prediction of all required coupons

%% set global variables
global test_strain test_stress % set the test strain and test stress as global variables

%% calculate the stress according to the test strain
simu_stress = monotonic_hardening(optimized_mono_mat_paras, test_strain);% compute simulation stress for a set of given material parameters

%% calculate the error
err_norm2 = norm((simu_stress(10:end) - test_stress(10:end))./test_stress(10:end)); % compute the difference

end
