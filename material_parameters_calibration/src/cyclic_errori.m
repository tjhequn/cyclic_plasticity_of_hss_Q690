function error = cyclic_errori(mat_paras,name_of_loading)

% This function is used to compute the difference of test stress and
% simulation stress.

%% generate the variable name
test_strain_name = [name_of_loading, '_test_strain']; % define the name of test strain
test_stress_name = [name_of_loading, '_test_stress']; % define the name of test stress
simu_stress_name = [name_of_loading, '_simu_stress']; % define the name of simulation stress

%% define global variable
eval(['global ', test_strain_name,' ', test_stress_name]) % set the test strain and test stress as global variables
global num_alpha ... % the number of kinematic hardening
       num_Y_iso ... % the number of isotropic hardening

%% calculate the stress according to the test strain
eval(['[', simu_stress_name, '] = cyclic_softening(mat_paras,num_alpha,num_Y_iso,', test_strain_name, ');']); % compute simulation stress for a set of given material parameters

%% calculate the error
error = eval(['norm(', simu_stress_name, ' - ', test_stress_name, ');']); % compute the difference

end
