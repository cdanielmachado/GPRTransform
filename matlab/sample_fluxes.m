%% Run sampling for GPR-transformed version of the ecoli core model

model = readCbModel('../models/ecoli_core_unfolded.xml');

changeCobraSolver('gurobi5', 'all');

points = 1e4;
steps = 1e2;

% run for wild-type

sol_wt = optimizeCbModel(model);
bio = strmatch('Biomass', model.rxns);
atpm = strmatch('ATPM', model.rxns);
model.ub(atpm) = model.lb(atpm);

model.lb(bio) = sol_wt.f * 0.9;
[reactions, sample] = sampler(model, points, steps);
save('../results/sampling/wt.mat', 'reactions', 'sample');

% run for succinate over-producer

succ = strmatch('EX_succ_e_f', model.rxns);
model.c([bio, succ]) = [0, 1];
model.lb(bio) = sol_wt.f*0.1;
sol_mut = optimizeCbModel(model);
model.lb(succ) = sol_mut.f*0.9;
[reactions, sample] = sampler(model, points, steps);
save('../results/sampling/succ.mat', 'reactions', 'sample');

% clean up

delete('ACHR_last_point.mat', 'ACHRerror.txt', 'modelRedTmp.mat', ...
    'sampleCbModelTmp.mat', 'samplingtempfiles_1.mat');



