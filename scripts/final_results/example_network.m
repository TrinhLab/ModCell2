%% Parse input and create prodnet
modcell_path            = fileparts(which('initModCell2.m'));
input_info.problem_path = fullfile(modcell_path,'problems','example_network');
prodnet                 = Prodnet(input_info); 
prodnet.save();

%% Design
prodnet.set_deletion_type('reactions');
de = MCdesign(prodnet);
% change default ga parameters:
de.ga_parameters.random_num_gen_seed = 1;
de.ga_parameters.use_parallel        = false;
de.ga_parameters.stall_generations   = 30;
de.ga_parameters.population_size     = 30;

%%
% wGCP objective and different parameters
design_parameters.objective     = 'wGCP';
design_parameters.max_deletions = 1;

% no module reactions
design_parameters.max_module = 0 .* ones(prodnet.n_prod,1);
de.solve_mop(design_parameters);

% one module reaction per network
design_parameters.max_module = 1 .* ones(prodnet.n_prod,1);
de.solve_mop(design_parameters);

% sGCP and different parameters
design_parameters.objective = 'sGCP';
design_parameters.max_deletions = 3;

% no module reactions
design_parameters.max_module = 0 .* ones(prodnet.n_prod,1);
de.solve_mop(design_parameters);

% one module reaction per network
design_parameters.max_module = 1 .* ones(prodnet.n_prod,1);
de.solve_mop(design_parameters);

%%
% NGP objective setup
design_parameters.objective = 'NGP';
prodnet.set_deletion_type('reactions','NGP'); % (This points to the same prodnet as prodnet in the "de" object.
design_parameters.max_deletions = 3;

% no module reactions
design_parameters.max_module = 0 .* ones(prodnet.n_prod,1);
de.solve_mop(design_parameters);

% one module reaction per network
design_parameters.max_module = 1 .* ones(prodnet.n_prod,1);
de.solve_mop(design_parameters);

%% Analyze results
%%
prodnet = load_prodnet('example_network');

%%
% Write solution reports
ra = ResAnalysis(prodnet,{'wGCP-1-0','wGCP-1-1','sGCP-3-0','sGCP-3-1'});
ra.write_to_xls('growth_obj_report')

ra = ResAnalysis(prodnet,{'NGP-3-0','NGP-3-1'});
ra.write_to_xls('non_growth_obj_report')

%%
% Plot production envelope matrix
ra = ResAnalysis(prodnet,{'wGCP-1-0','wGCP-1-0','wGCP-1-1',...
    'sGCP-3-0','sGCP-3-0','sGCP-3-0','sGCP-3-1',...
    'NGP-3-0','NGP-3-0','NGP-3-1'});

ra.plot_yield_vs_growth([1,2,1,...
    1,2,3,1,...
    1,2,1])
