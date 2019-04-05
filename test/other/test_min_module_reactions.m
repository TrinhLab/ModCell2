function test_min_module_reactions(ra, sol_ind, des_ind)
% Test for minimize designs
%
% Args:
%   - ra (ResAnalysis) : Instance of result analysis that has been
%   subject to min_module_reactions. 
%   - sol_ind (integer): Solution index
%   - des_ind (integer): Design index
%   
% Notes:
%   -module minimization might lead to small discrepancy due to
% objective tolerance

if ~exist('OBJ_TOL', 'var')
	OBJ_TOL = 1e-5;
end
original_obj = ra.solutions(sol_ind).design_objectives(des_ind,:)';

y = ra.solutions(sol_ind).design_deletions(des_ind,:);
Z = ra.solutions(sol_ind).design_modules(des_ind).Z;
ra.prodnet.set_module_and_deleted_variables(Z,y)
cur_obj = ra.prodnet.calc_design_objectives(ra.solutions(sol_ind).design_parameters.objective);

assert (all( abs(original_obj - cur_obj) <= OBJ_TOL))
end