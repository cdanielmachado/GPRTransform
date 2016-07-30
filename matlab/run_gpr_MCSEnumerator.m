%% run MCSEnumerator on a GPR-transformed model

cbmodel = readCbModel('../models/iAF1260_MCSEnum_unfolded.xml');

max_dels = 8;
substrate = 'EX_glc_e_b';
product = 'EX_etoh_e_f';
min_yield = 1.4;
min_growth = 0.001;
glc_uptake = 10; 
cbmodel.ub(strcmp(cbmodel.rxns, substrate)) = glc_uptake;

notknockable=find(not(strncmp('u_', cbmodel.rxns, 2)))';

cmcs = MCSEnumerator_wrapper(cbmodel, max_dels, substrate, product, min_yield, min_growth, notknockable);

save('../results/mcs/mcs_iAF1260_gpr_MCSEnum_max8.mat', 'cmcs') 



