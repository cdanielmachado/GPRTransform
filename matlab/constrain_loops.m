function model_new = constrain_loops(model)
% Utility function to constrain futile cycles

objPercentage = 0;

k = length(model.genes) - 1;
[m, n] = size(model.S);
model2.S = model.S(1:m-k, 1:n-k);
model2.rev = model.rev(1:n-k);
model2.lb = model.lb(1:n-k);
model2.ub = model.ub(1:n-k);
model2.c = model.c(1:n-k);
model2.rxns = model.rxns(1:n-k);
model2.mets = model.mets(1:m-k);

[minFlux, maxFlux] = fluxVariability(model2, objPercentage, 'max', model2.rxns, false, false);

model_new = model;

model_new.lb(1:n-k) = minFlux;
model_new.ub(1:n-k) = maxFlux;

end


