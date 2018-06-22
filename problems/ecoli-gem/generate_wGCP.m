function generate_wGCP(max_deletions, max_module, MAX_TIME)
modcell_path = fileparts(which('initModCell2.m'));
%%
diary(['wGCP_',num2str(max_deletions),'.log'])
lin = load(fullfile(modcell_path,'problems','ecoli-gem', 'prodnet-known-l.mat'));
prodnet = lin.prodnet;
prodnet.problem_path = fullfile(modcell_path,'problems',prodnet.problem_name);

prodnet.set_deletion_type('reactions');
de = MCdesign(prodnet);

% change default ga parameters:
de.ga_parameters.stall_generations   = 500;
de.ga_parameters.population_size     = 400;

% start point
start_id = ['wGCP-', num2str(max_deletions-1), '-', num2str(max_module)];
de.set_start_point(start_id,prodnet.problem_path);

% design parameters
design_parameters.objective     = 'wGCP';
design_parameters.max_deletions = max_deletions;
design_parameters.max_module    = max_module .* ones(prodnet.n_prod,1);
max_gen = 10000;

% solve
de.solve_mop(design_parameters, [], max_gen, MAX_TIME);
diary off
end
