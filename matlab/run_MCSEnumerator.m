%% Run MCSEnumerator on iAF1260 (reproduces results from the original paper)

cbmodel = readCbModel('../models/iAF1260_MCSEnum.xml');

load('notknockable.mat') %these are defined in the original paper

max_dels = 7;
substrate = 'EX_glc(e)';
product = 'EX_etoh(e)';
oxygen = 'EX_o2(e)';
min_yield = -1.4;
min_growth = 0.001;
glc_uptake = 10; 
cbmodel.lb(strcmp(cbmodel.rxns, oxygen)) = 0;
cbmodel.lb(strcmp(cbmodel.rxns, substrate)) = -glc_uptake;

cmcs = MCSEnumerator_wrapper(cbmodel, max_dels, substrate, product, min_yield, min_growth, notknockable);

save('../results/mcs/mcs_iAF1260_MCSEnum_max7.mat', 'cmcs') 


