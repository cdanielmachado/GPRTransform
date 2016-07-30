function [reactions, sample] = sampler(model, points, steps)
% Simple wrapper for the cobra ACHR sampler

model = constrain_loops(model);

method = 'ACHR';
output = 'samplingtempfiles';
options.nWarmupPoints = 1e4;
options.nFiles = 1;
options.nPointsPerFile = points;
options.nStepsPerPoint = steps;
options.nPointsReturned = points;
options.nFilesSkipped = 0;
options.removeLoopsFlag = false;
options.removeLoopSamplesFlag = false;

[redmodel,sample] = sampleCbModel(model,output,method,options);
reactions = redmodel.rxns;

end

