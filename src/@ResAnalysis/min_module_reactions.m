function min_module_reactions(obj, OBJ_TOL)
% Removes unnecessary module reactions from soluitons in the ResAnalysis object.
%
% Args:
%	OBJ_TOL: Minimum value for an objective to be considered differnt (default 1e-5)
%
%
% Notes
% -----
% - parallelize for each solution?
% - create a wraper around fprintf that only prints if global verbose variable is active (similar to decorator)"
% - create a core function that is as indepent as possible of the ResAnalysis data structure


if ~exist('OBJ_TOL', 'var')
	OBJ_TOL = 1e-5;
end

fprintf('module minimization..\n')
obj.prodnet.reset_wild_type_state();
obj.prodnet.set_deletion_type('reactions');
for sol_ind = 1:length(obj.solutions)
	if isempty(obj.solutions(sol_ind).design_modules)
		continue;
	end
	dp = obj.solutions(sol_ind).design_parameters;
	fprintf('\t %s-%d-%d (%d)\n',dp.objective, dp.max_deletions,max(dp.max_module), sol_ind)

	fprintf('\t\t design  \t product \t old_modules  \t --> new_modules\n')
	for des_ind = 1:length(obj.solutions(sol_ind).design_modules)
		y = obj.solutions(sol_ind).design_deletions(des_ind,:);
		Z = obj.solutions(sol_ind).design_modules(des_ind).Z;
		best_obj = obj.solutions(sol_ind).design_objectives(des_ind,:)';

		newZ = false(size(Z,1), size(Z,2));
		for prod_ind = 1:obj.prodnet.n_prod
			if sum(Z(prod_ind,:)) == 0 % Ignore designs not using modules for this prod
				continue;
			end

			% Try combinations of increasing cardinality until objective matches, if no match then leave original modules
			mod_idx = find(Z(prod_ind, :));
			found_sol = false;
			for card = 0:length(mod_idx) - 1

				if found_sol
					break;
				end

				combs = combnk(mod_idx, card);
				for comb_ind = 1:size(combs, 1)
					tempZ = false(size(Z,1), size(Z,2));
                    if card > 0
                        tempZ(prod_ind, combs(comb_ind)) = true;
                    end
					obj.prodnet.set_module_and_deleted_variables(tempZ,y)
					cur_obj = obj.prodnet.calc_design_objectives(obj.solutions(sol_ind).design_parameters.objective, prod_ind);
					if abs(cur_obj(prod_ind) - best_obj(prod_ind)) < OBJ_TOL
						found_sol = true;
					break;
					end
				end
			end
			if found_sol
				newZ(prod_ind,:) = tempZ(prod_ind,:);
				old_mod = join(obj.prodnet.parent_model.rxns(obj.prodnet.cand_ind(Z(prod_ind,:))), ',');
				new_mod = join(obj.prodnet.parent_model.rxns(obj.prodnet.cand_ind(newZ(prod_ind,:))), ',');
				fprintf('\t\t %d \t %s \t %s \t --> %s\n', des_ind, obj.prodnet.prod_id{prod_ind}, old_mod{1}, new_mod{1})
			else
				newZ(prod_ind,:) = Z(prod_ind,:);
			end
		end
		% Overwrite
		obj.solutions(sol_ind).design_modules(des_ind).Z = newZ;
	end
end
fprintf('..module minimization done\n')
end
