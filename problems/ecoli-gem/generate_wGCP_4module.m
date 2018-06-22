clear; clc
modcell_path = fileparts(which('initModCell2.m'));
diary wGCP-4-module.log
%% 4-1,2,3
for module_ind = 1:3
load(fullfile(modcell_path,'problems','ecoli-gem', 'prodnet-known-l.mat'))
prodnet.problem_path = fullfile(modcell_path,'problems',prodnet.problem_name);
prodnet.LP_SOLVER = 'glpk';
prodnet.set_deletion_type('reactions');

prodnet.set_deletion_type('reactions');
de = MCdesign(prodnet);

% change default ga parameters:
de.ga_parameters.stall_generations   = 500;
de.ga_parameters.population_size     = 400;

de.set_start_point(['wGCP-4-',num2str(module_ind-1)],prodnet.problem_path)
% design parameters
design_parameters.objective     = 'wGCP';
design_parameters.max_deletions = 4;
design_parameters.max_module    = module_ind .* ones(prodnet.n_prod,1);
max_gen = 10000;
max_time = 15*60; % 15 hours

% solve
de.solve_mop(design_parameters, [], 10000, max_time);
end
diary off
