%% compute EFMs for simple glycolysis model using CellNetAnalyser

cbmodel = readCbModel('../models/glycolysis.xml');
cnap = CNAcobra2cna(cbmodel);
cnap.reacMin(cnap.reacMin == -1000) = -inf;
cnap.reacMax(cnap.reacMax == 1000) = inf;

[efms,rev,idx,ray] = CNAcomputeEFM(cnap);

save('../results/EFMs/EFMs_glycolysis.mat', 'efms');

cbmodel = readCbModel('../models/glycolysis_unfolded.xml');
cnap = CNAcobra2cna(cbmodel);
cnap.reacMin(cnap.reacMin == -1000) = -inf;
cnap.reacMax(cnap.reacMax == 1000) = inf;

[efms,rev,idx,ray] = CNAcomputeEFM(cnap);

save('../results/EFMs/EFMs_glycolysis_unfolded.mat', 'efms');