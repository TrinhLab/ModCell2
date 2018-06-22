clear; clc;
modcell_path = fileparts(which('initModCell2.m'));
diary create_prodnet.log

%% all b
input_info.problem_path = fullfile(modcell_path,'problems','ecoli-gem'); 
input_info.parent_model_path = fullfile(fullfile(modcell_path,'problems','ecoli-gem-set48'), 'parent', 'parent-all-b-48.mat');
prodnet = Prodnet(input_info); 
prodnet.save('prodnet-all-b');

diary off