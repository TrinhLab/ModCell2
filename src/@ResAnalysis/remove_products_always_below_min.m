function remove_products_always_below_min(obj, min_obj)
% Removes products that do not attain a minimum objective in any solution
%
% Args:
%   min_obj (float): At least one design with objective value greater or equal than min_obj should exist for each product that is maintained in the ResAnalysis structure.

above_min = false(1,length(obj.prodnet.prod_id));
for sol_ind = 1:length(obj.solutions)
    above_min = above_min | any(obj.solutions(sol_ind).design_objectives >= min_obj,1);
end
prod_to_remove = obj.prodnet.prod_id(~above_min);
fprintf('The following products were removed under the minimum objective %2.2f:\n', min_obj) % remove_products will print them
obj.remove_products(prod_to_remove);
end
