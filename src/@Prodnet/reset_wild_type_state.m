function reset_wild_type_state(obj)
% Sets production networks to their wild type state.
%

for i =1: obj.n_prod
    obj.model_array(i).c(:) = 0;
    obj.model_array(i).lb = obj.model_array(i).original_lb;
    obj.model_array(i).ub = obj.model_array(i).original_ub;
    obj.model_array(i).module_rxn_ind = [];
    obj.model_array(i).module_gene_ind = [];
end