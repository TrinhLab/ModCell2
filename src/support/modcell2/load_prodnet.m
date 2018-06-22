function prodnet = load_prodnet(problem_name)
% Load production network ensuring that the problem_path field is correct
% 
% Args:
%    problem_name(string)

modcell_path = fileparts(which('initModCell2.m'));
prodnet_path = fullfile(modcell_path,'problems',problem_name,'prodnet.mat');

loadin  = load(prodnet_path);
prodnet = loadin.prodnet;

% Set appropriate path
prodnet.problem_path = fullfile(modcell_path,'problems',prodnet.problem_name);

end