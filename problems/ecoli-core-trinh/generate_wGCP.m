function generate_wGCP(max_deletions, max_module, MAX_TIME)
modcell_path = fileparts(which('initModCell2.m'));

problem_id = 'ecoli-core-trinh';
%%
diary(['wGCP_', num2str(max_deletions), '_', num2str(max_module), '.log'])
lin = load(fullfile(modcell_path,'problems',problem_id, 'prodnet.mat'));
prodnet = lin.prodnet;
prodnet.problem_path = fullfile(modcell_path,'problems',prodnet.problem_name);

prodnet.LP_SOLVER = 'cplex';

prodnet.set_deletion_type('reactions');
de = MCdesign(prodnet);

% change default ga parameters:
de.ga_parameters.stall_generations   = 300;
de.ga_parameters.population_size     = 200;

% start point
if max_deletions > 1
    start_id = ['wGCP-', num2str(max_deletions-1), '-', num2str(max_module)];
try
    de.set_start_point(start_id,prodnet.problem_path);
catch
    fprintf('Start point not found: %s\n',start_id);
end
% design parameters
design_parameters.objective     = 'wGCP';
design_parameters.max_deletions = max_deletions;
design_parameters.max_module    = max_module .* ones(prodnet.n_prod,1);
max_gen = 10000;

% solve
de.solve_mop(design_parameters, [], max_gen, MAX_TIME);

poolobj = gcp('nocreate');
delete(poolobj);
diary off
end
