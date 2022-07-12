function plot_comparison(mat_paras,num_alpha,num_Y_iso,coupon_name,...
                         test_strain,test_stress,export_data_save_switch,export_data_save_path)

% This function is used to plot comparison between test results and simulation results
    
%% calculate simulation results
simu_strain = test_strain;
[simu_stress,alpha,Y_iso,estrain,pstrain,peeq,psi,phi,E_modulus,num_iteration] = ... 
    cyclic_softening(mat_paras,num_alpha,num_Y_iso,simu_strain); % calculate the simulation stress for the given test strain
relative_err = (simu_stress-test_stress)/max(test_stress)*100; % calculate the relative error

%% plot figures for test validation
% show the complete difference
figure_name = [coupon_name, '-full'];
fig_comparison_full = figure('NumberTitle', 'off', 'Name', figure_name); %define the figure name
set(gcf, 'position', [200,100,1360,960]) % define the figure size

subplot(2,2,1) % comparison of stress strain response
plot(test_strain*100,test_stress,'linewidth',1)
hold on
plot(simu_strain*100,simu_stress, '--r','linewidth',1)
legend('test','simulation', 'Fontsize', 12,'Location','NorthWest')
title('Hysteretic loops', 'Fontsize', 14)
xlabel('Strain (%)') %x-axis
ylabel('Stress (MPa)') %y-axis

subplot(2,2,2) % show the relative error
histogram(relative_err,'Normalization','probability')
title('Relative error with maximum stress', 'Fontsize', 14)
xlabel('Err (%)') %x-axis
ylabel('P') %y-axis

subplot(2,2,[3,4]) % comparison of peeq stress comparison
plot(peeq,test_stress,'Linewidth',1)
hold on
plot(peeq,simu_stress, '--r','Linewidth',1)
legend('test','simulation', 'Fontsize', 12,'Location','NorthEast')
title('peeq-stress', 'Fontsize', 14)
xlabel('PEEQ') %x-axis
ylabel('Stress (MPa)') %y-axis

% separated cycles
figure_name = [coupon_name, '-separated'];
fig_comparison_separated = figure('NumberTitle', 'off', 'Name', figure_name); % define the figure name
set(gcf, 'position', [200,100,1360,960]) % define the figure size

num_x = 4; % set the number of subplot in horizontal direction
num_y = 4; % set the number of subplot in c direction
total_num = num_x*num_y; % compute the total number of subplot
num_data_each = floor(numel(test_strain)/total_num); % equally divide the test data

for j = 1:total_num % loop over all subplots
    start_id = (j-1)*num_data_each + 1;
    end_id   = min(j*num_data_each,numel(test_strain));
    temp_test_strain = test_strain(start_id:end_id);
    temp_test_stress = test_stress(start_id:end_id);
    temp_simu_stress = simu_stress(start_id:end_id);

    subplot(num_x,num_y,j)
    plot(temp_test_strain*100,temp_test_stress,'linewidth',1)
    hold on
    plot(temp_test_strain*100,temp_simu_stress, '--r','linewidth',1)
    legend('test','simulation', 'Fontsize', 8,'Location','NorthWest')
    title(['part-', num2str(j)], 'Fontsize', 10)
    xlabel('Strain (%)') %x-axis
    ylabel('Stress (MPa)') %y-axis

end

% export data and save it
if strcmp(export_data_save_switch,'on') || strcmp(export_data_save_switch,'ON')
    column_names = {   'test_strain','test_stress',  'simu_strain','simu_stress','alpha','Y_iso','peeq','relative_err'};
    tmp_table    = table(100*test_strain,  test_stress,100*simu_strain,  simu_stress,  alpha,  Y_iso,  peeq,  relative_err,'VariableNames',column_names);
    writetable(tmp_table,[export_data_save_path, coupon_name,'_simulation_results.csv']);
    
    file_name_full = [export_data_save_path, coupon_name, '-full'];
    saveas(fig_comparison_full, file_name_full, 'png') % save the figure as .png file
    
    file_name_separated = [export_data_save_path, coupon_name, '-separated'];
    saveas(fig_comparison_separated, file_name_separated, 'png') % save the figure as .png file
else
    disp('No figures are saved! If you want to save figures, please set export_data_save_switch = ''on'' or ''ON''!')
end

end
