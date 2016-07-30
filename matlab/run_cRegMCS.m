%% Run cRegMCS on a GPR-transformed ecoli core model

cbmodel = readCbModel('../models/ecoli_core_unfolded.xml');

changeCobraSolver('gurobi5');

substrate = 'EX_glc_e_b';
product = 'EX_etoh_e_f';
biomass = 'Biomass_Ecoli_core_w_GAM';
oxygen = 'EX_o2_e_b';

substrate_idx = find(strcmp(cbmodel.rxns, substrate));
product_idx = find(strcmp(cbmodel.rxns, product));
biomass_idx = find(strcmp(cbmodel.rxns, biomass));
oxygen_idx = find(strcmp(cbmodel.rxns, oxygen));

max_dels = 3;
min_growth = 0.001;
min_yield = 1.4;

cbmodel.ub(oxygen_idx) = 0;

notknockable=find(not(strncmp('u_', cbmodel.rxns, 2)))';

regulated = find(strncmp('u_', cbmodel.rxns, 2))';

coupled = {'u_b0733', 'u_b0734', 'u_b0809', 'u_b0810', 'u_b0811', 'u_b0978', 'u_b0979'};

coupled_idx = zeros(1, length(coupled));
for i = 1:length(coupled)
    coupled_idx(i) = find(strcmp(cbmodel.rxns, coupled{i}));
end
regulated = setdiff(regulated, coupled_idx);

[vmin, vmax] = fluxVariability(cbmodel,  0, 'max', cbmodel.rxns(regulated));
regulation.reg_ind = regulated(vmax < 500); 
regulation.numregsteps=3; 
regulation.reg_down_up=[true(1,numel(regulation.reg_ind));
                        true(1,numel(regulation.reg_ind))];
regulation.reg_bounds=[];

cnap = CNAcobra2cna(cbmodel);
cnap.reacMax(cnap.reacMax >= 1000) = inf;
cnap.reacMin(cnap.reacMin <= -1000) = -inf;


n = length(cbmodel.rxns);
T = zeros(1, n);
T(1, product_idx) = 1; %#ok<*FNDSB>
T(1, substrate_idx) = -min_yield;
t = 0;

biomass = find(cbmodel.c);
D= zeros(1, n);
D(1, biomass)= -1;
d= -min_growth;

outputfile = 'regmcs_etoh.txt';

cmcs = CNAregMCSEnumerator(cnap, T, t, D, d, notknockable, inf, max_dels, outputfile,0,[],regulation);



