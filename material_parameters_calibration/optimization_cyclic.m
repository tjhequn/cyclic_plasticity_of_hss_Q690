% This program is used to calibrate material parameters for cyclic loading
% via optimization method.

close all
clear
clc

tic % mark the start time

%% temporarily add the path of source code
addpath(genpath('./src'));

%% set global variables
global optimization_data_name_set ... % coupon name list for optimization
       num_alpha ... % the number of kinematic hardening
       num_Y_iso ... % the number of isotropic hardening
       elastic_paras ... % elastic parameters, i.e., elastic modulus
       monotonic_paras ... % material parameters for monotonic loading
       elastic_modulus_paras ... % parameters for elastic modulus evolution

%% set the dataset for optimization
test_data_file_path  = '.\test_data\'; % define the test data file path
all_coupon_name_list = {'CTS_01','CTS_02','CTS_03','CTS_04','CTS_05','CTS_06','CTS_07',...
                        'CTS_08','CTS_09','CTS_10','CTS_11','CTS_12','CTS_13'}; % define all coupon names
require_coupon_index = [6,7]; % define the coupon required for optimization
optimization_data_folder = '.\test_data\'; % define the folder for saving the test data
optimization_data_name_set = all_coupon_name_list(require_coupon_index);

%% set lower and upper bound for optimization
m_phi_lower_upper        = [0.0,100.0];
phi_sat_lower_upper      = [350,550];
m_alpha_list_lower_upper = [1000,10000; 100,1000; 10,100];
omega_list_lower_upper   = [0.0,1.0; 0.0,1.0; 0.0,1.0];
m_y_list_lower_upper     = [1000,5000; 1.49,1.51; 0.049,0.051];
Q_list_lower_upper       = [-391,-380; -56,-55; -101,-100];

num_alpha = numel(m_alpha_list_lower_upper(:,1));
num_Y_iso = numel(m_y_list_lower_upper(:,1));

%% basic material parameters, including elastic constants and monotonic parameters
% mataPara is a column vector that stores material parameters.
% elastic parameters
E_mod0        = 204255; % elastic modulus
elastic_paras = E_mod0; % generate elastic parameter vector

% monotonic plastic parameters
sigma_y         = 799.8; % yield stress
strain_sh       = 0.010; % length of yield plateau
sigma_sat       = 247.8; % nonlinear hardening
m_N             = 7.67;  % nonlinear hardening modulus
m_L             = 249.8; % linear hardening modulus
monotonic_paras = [sigma_y, strain_sh, sigma_sat, m_N, m_L]; % generate monotonic prameter vector

% elastic modulus parameters
E_mod_sat             = 178308;
xi_E                  = 95.5;
elastic_modulus_paras = [E_mod_sat, xi_E];

%% read test data
for i  = 1:numel(optimization_data_name_set)
    
    num_interval = 4; % define the data points to be skipped
    name_of_loading = optimization_data_name_set{i};
    load([optimization_data_folder, name_of_loading, '.mat'])

    test_strain_name = [name_of_loading, '_test_strain']; % define the name of test strain
    test_stress_name = [name_of_loading, '_test_stress']; % define the name of test stress
    simu_stress_name = [name_of_loading, '_simu_stress']; % define the name of simulation stress

    eval(['global ', test_strain_name,' ', test_stress_name]) % set the test strain and test stress as global variables
    eval([test_strain_name, ' = test_data(1:num_interval:end,2)/100;']); % read test strain
    eval([test_stress_name, ' = test_data(1:num_interval:end,3);']); % read test stress
    
    eval(['clear ', name_of_loading]); % clear the name of test data
    clear name_of_loading simu_stress_name test_strain_name test_stress_name % clear the names defined above
end

clear i num_interval % clear unnecessary variables

%% intialize parameters for optimization
% generate the linear inequality constraint
A = []; % linear inequality constraint
b = []; % linear inequality constraint

% generate the linear equality constraint
Aeq = zeros(1,2+2*num_alpha+2*num_Y_iso); % linear equality constraint
Aeq(3+num_alpha:2+2*num_alpha) = ones(1,num_alpha);
beq = [1]; % linear equality constraint

% generate the lower bound and upper bound for parameters
lb = [m_phi_lower_upper(:,1); phi_sat_lower_upper(:,1); m_alpha_list_lower_upper(:,1); ...
      omega_list_lower_upper(:,1); m_y_list_lower_upper(:,1); Q_list_lower_upper(:,1)]'; % lower bound
ub = [m_phi_lower_upper(:,2); phi_sat_lower_upper(:,2); m_alpha_list_lower_upper(:,2); ...
      omega_list_lower_upper(:,2); m_y_list_lower_upper(:,2); Q_list_lower_upper(:,2)]';  % upper bound
nonlcon = []; % nonlinear constraint

% set the initial value for optimization
x0 = 0.5*lb + 0.5*ub; % intial point for optimization
x0(3+num_alpha:2+2*num_alpha) = ones(1,num_alpha)/num_alpha; % initialize omega_i

%% set conditions for optimization control such as optimization algorithm, stopping criteria
options = optimoptions('fmincon','Display','iter','Algorithm','sqp'); % set the optimization algorithm
optim_paras = fmincon(@cyclic_error_calculation,x0,A,b,Aeq,beq,lb,ub,[],options); % start optimization, the objective function is cyclic_error_calculation

%% generate material parameter vector and export optimized material parameters
% generate the material parameter vector
optimized_mat_paras = [elastic_paras, monotonic_paras, optim_paras, elastic_modulus_paras];

% export the optimized parameters into a table
row_names = {'E';'sigma_y';'strain_sh';'sigma_sat';'m_N';'m_L';'m_phi';'phi_sat'};
m_alpha_name_list = repmat({'m_alpha'},[num_alpha,1]);
omega_name_list   = repmat({'omega'},[num_alpha,1]);
m_y_name_list     = repmat({'m_y'},[num_Y_iso,1]);
Q_i_name_list     = repmat({'Q_i'},[num_Y_iso,1]);

for i = 1:num_alpha
    m_alpha_name_list{i} = ['m_alpha_',num2str(i)];
    omega_name_list{i}   = ['omega_',num2str(i)];
end

for i = 1:num_Y_iso
    m_y_name_list{i} = ['m_y_',num2str(i)];
    Q_i_name_list{i} = ['Q_i_',num2str(i)];
end

parameter_names   = [row_names;m_alpha_name_list;omega_name_list;m_y_name_list;Q_i_name_list;{'E_sat';'xi_E'}];
optimized_values  = optimized_mat_paras';
optimization_data = table(parameter_names,optimized_values);
writetable(optimization_data,'optimized_material_parameters.csv');

%% plot simulation results and save figures
for i  = 1:numel(require_coupon_index)
    % read test data
    required_coupon_number = require_coupon_index(i);
    coupon_name = all_coupon_name_list{required_coupon_number}; % get the coupon name from the coupon name list
    test_data_file_name = [test_data_file_path,coupon_name]; % define the file of test data
    [test_strain,test_stress] = read_test_data(test_data_file_name); % read the test data and get the test strain and test stress
    
	% plot figures and do not export any figures and data
    export_data_save_switch = 'off';
    export_data_save_path   = 'null';
    plot_comparison(optimized_mat_paras,num_alpha,num_Y_iso,coupon_name,test_strain,test_stress,...
                    export_data_save_switch,export_data_save_path);
end

toc % mark the stop time
