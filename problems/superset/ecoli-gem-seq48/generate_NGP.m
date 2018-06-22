function generate_NGP(max_deletions, max_module, MAX_TIME, start_id)
if ~exist('start_id')
start_id = ['NGP-', num2str(max_deletions-1), '-', num2str(max_module)];
end
modcell_path = fileparts(which('initModCell2.m'));

%%
diary(['NGP_', num2str(max_deletions), '_', num2str(max_module), '.log'])
lin = load(fullfile(modcell_path,'problems','ecoli-gem-seq48', 'prodnet-all-b.mat'));
prodnet = lin.prodnet;
prodnet.problem_path = fullfile(modcell_path,'problems','ecoli-gem-seq48');

prodnet.LP_SOLVER = 'linprog';

prodnet.set_deletion_type('reactions', 'NGP');
de = MCdesign(prodnet);

% change default ga parameters:
de.ga_parameters.stall_generations   = 500;
de.ga_parameters.population_size     = 400;

% start point
try    
    de.set_start_point(start_id,prodnet.problem_path);
catch
fprintf('Start point not available\n')
end
% design parameters
design_parameters.objective     = 'NGP';
design_parameters.max_deletions = max_deletions;
design_parameters.max_module    = max_module .* ones(prodnet.n_prod,1);
max_gen = 10000;

% solve
de.solve_mop(design_parameters, [], max_gen, MAX_TIME);
poolobj = gcp('nocreate');
delete(poolobj);
diary off
end
