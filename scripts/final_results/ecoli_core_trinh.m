%% Generate e.coli core fixed module reactions:

%% 
% collect input data
input_info.problem_path = 'C:\Users\sergio\Dropbox\modcell2-software\problems\ecoli-core-trinh';
prodnet                 = Prodnet(input_info); 
prodnet.save();

%% Design modular cells
prodnet.set_deletion_type('reactions');
de = MCdesign(prodnet);
de.ga_parameters.use_parallel           = true;

% change default ga parameters:
de.ga_parameters.stall_generations      = 500;
de.ga_parameters.population_size        = 200;
de.ga_parameters.random_num_gen_seed    = 999; % Results generated in parallel MAY NOT be reproducible with this method

design_parameters.objective     = 'sGCP';
design_parameters.max_deletions = 5;
design_parameters.max_module    = 1 .* ones(prodnet.n_prod,1);

de.solve_mop(design_parameters);

%enumerate alternative solutions:
% targets in each case (total number known from MIP solution)

% Other parameters:

design_parameters.max_deletions = 6;
design_parameters.max_module    = 1 .* ones(prodnet.n_prod,1);
de.solve_mop(design_parameters);

design_parameters.max_deletions = 8;
design_parameters.max_module    = 2 .* ones(prodnet.n_prod,1);
de.solve_mop(design_parameters);

%% plot production spaces for the solutiosn matching those obtained with the MIP method
ra = Result_analysis(prodnet, {'sGCP-5-1','sGCP-6-1','sGCP-8-2'});
ra.plot_yield_vs_growth2([8,10,9],'convHullLineWidthMut',3,'yticks_values',[0,0.6]);

%% write output to excel:
ra.write_to_xls();

%% ignore
load('C:\Users\Sergio Garcia\Dropbox\modcell2-software\problems\ecoli-core-trinh\prodnet.mat')
set_prodnet(prodnet,1)

