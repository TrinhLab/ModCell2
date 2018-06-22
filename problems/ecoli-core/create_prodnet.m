tic
clear; clc;
modcell_path = fileparts(which('initModCell2.m'));
input_info.problem_path =fullfile(modcell_path,'problems','ecoli-core'); 
prodnet = Prodnet(input_info); 
prodnet.problem_path = []; % Avoid forgetting to set this
prodnet.save('prodnet');
toc