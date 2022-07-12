function [test_strain,test_stress] = read_test_data(test_data_file_name)
% This function is used to read the strain and stress histroy from test
% data, it is just a sample code. Users can rewrite it according to their
% demand.
% I have made a pre-processing of my test data. All of them have been saved
% in .mat format. The first column is the loading cycle, the second column
% is the strain with the unit of percentage, the third column is the stress
% with the unit MPa. 

load([test_data_file_name, '.mat'])
test_strain = test_data(:,2)/100; % for the strain in precentage
test_stress = test_data(:,3);

end