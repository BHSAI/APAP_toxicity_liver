function [timbr_network_demand, timbr_parsimonious_solution] = ...
    timbr(model, production_rxn, production_required, timbr_weights)

if (strcmp(production_rxn,''))
    timbr_model = model;
else
    timbr_model = changeRxnBounds(model,production_rxn,production_required,'l');
    timbr_model = changeRxnBounds(timbr_model,production_rxn,1000000,'u');
end

reverse_rxn = regexprep(production_rxn,'_f','_r');
if (any(strcmp(reverse_rxn,timbr_model.rxns)) && ...
        ~any(strcmp(reverse_rxn,production_rxn)) && ...
        any(strcmp(production_rxn,timbr_model.rxns)))
     if (strcmp(reverse_rxn, 'RCR30327_r')) % glycerol exchange reversible =RCR30327_r
         timbr_model = changeRxnBounds(timbr_model,reverse_rxn,0,'l');
     else
    timbr_model = changeRxnBounds(timbr_model,reverse_rxn,0,'l');
    timbr_model = changeRxnBounds(timbr_model,reverse_rxn,0,'u');
     end
end

if (istable(timbr_weights))
    timbr_weights = table2array(timbr_weights);
end

timbr_model.c = timbr_weights;
timbr_parsimonious_solution = optimizeCbModel(timbr_model,'min');
timbr_network_demand = timbr_parsimonious_solution.f;

end
