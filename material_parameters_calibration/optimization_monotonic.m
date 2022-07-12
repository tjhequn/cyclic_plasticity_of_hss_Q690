% This program is used to calibrate material parameters for cyclic loading
% via optimization method.

close all
clear
clc

tic % mark the start time

%% temporarily add the path of source code
addpath(genpath('./src'));

%% set global variables
global coupon_name % coupon name list for optimization

%% set the dataset for optimization
test_data_file_path  = '.\test_data\'; % define the test data file path
coupon_name          = 'MTS_01'; % define all coupon names
optimization_data_folder = '.\test_data\'; % define the folder for saving the test data

%% set lower and upper bound for optimization
E_mod_lower_upper     = [204225,204225];
sigma_y_lower_upper   = [799.8,800];
strain_sh_lower_upper = [0.005,0.01];
sigma_sat_lower_upper = [0.0,500.0];
m_N_lower_upper       = [0,50];
m_L_lower_upper       = [0,500];

%% read test data  
num_interval = 1; % define the data points to be skipped
load([optimization_data_folder, coupon_name, '.mat'])

global test_strain test_stress % set the test strain and test stress as global variables
test_strain = test_data(1:num_interval:end,1); % read test strain
test_stress = test_data(1:num_interval:end,2); % read test stress

eval(['clear ', coupon_name]); % clear the name of test data
clear name_of_loading simu_stress_name test_strain_name test_stress_name % clear the names defined above

clear num_interval % clear unnecessary variables

%% intialize parameters for optimization
% generate the linear inequality constraint
A = []; % linear inequality constraint
b = []; % linear inequality constraint

% generate the linear equality constraint
Aeq = []; % linear equality constraint
beq = []; % linear equality constraint

% generate the lower bound and upper bound for parameters
lb = [E_mod_lower_upper(:,1); sigma_y_lower_upper(:,1); strain_sh_lower_upper(:,1); sigma_sat_lower_upper(:,1); ...
      m_N_lower_upper(:,1); m_L_lower_upper(:,1)]'; % lower bound
ub = [E_mod_lower_upper(:,2); sigma_y_lower_upper(:,1); strain_sh_lower_upper(:,2); sigma_sat_lower_upper(:,2); ...
      m_N_lower_upper(:,2); m_L_lower_upper(:,2)]';  % upper bound
nonlcon = []; % nonlinear constraint

% set the initial value for optimization
x0 = 0.5*lb + 0.5*ub; % intial point for optimization

%% set conditions for optimization control such as optimization algorithm, stopping criteria
options = optimoptions('fmincon','Display','iter','Algorithm','sqp'); % set the optimization algorithm
options.MaxFunctionEvaluations = 3000;
mono_optim_paras = fmincon(@monotonic_error_calculation,x0,A,b,Aeq,beq,lb,ub,[],options); % start optimization, the objective function is error_calculation

% export the optimized parameters into a table
parameter_names = {'E';'sigma_y';'strain_sh';'sigma_sat';'m_N';'m_L'};
optimized_values  = mono_optim_paras';
optimization_data = table(parameter_names,optimized_values);
writetable(optimization_data,'optimized_monotonic_material_parameters.csv');

%% plot simulation results and save figures
[simu_stress,R_iso,estrain,pstrain,peeq,num_iteration] = monotonic_hardening(mono_optim_paras,test_strain);
relative_err = (simu_stress-test_stress)/max(test_stress)*100;

%% plot figures
figure()
% set(gcf, 'position', [200,100,660,510])
plot(test_strain,test_stress)
hold on
plot(test_strain,simu_stress,'--r')
ylim([0,1000])
xlabel('strain')
ylabel('stress (MPa)')
legend('test','simulation','Location','southeast')

toc % mark the stop time
