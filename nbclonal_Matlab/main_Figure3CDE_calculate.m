function void = main_Figure3CDE_calculate(void)

clear all; close all; clc;

params.R0 = 1.6;
params.mulist = 0.01:0.01:2;
params.N0list = 1:8;

params.maxFinalSize = 120;
params.maxGeneratedMutantLineages = 200;
params.maxEstablishedMutantLineages = 200;
params.maxClonal = 10;

outfile3CDE = 'figure3CDE_data';

cntr = 1;
cntrN = 1;
for N0 = params.N0list
    params.N = N0;
    cntrmu = 1;
    for mu = params.mulist
        [N0 mu]
        params.mu = mu;
        results.muNvals(cntr,1:2) = [N0 mu];
        
        [n_clonals, rho] = calculateClonalPMF(params, 0);
        results.rhovals(cntr).n_clonals = n_clonals;
        results.rhovals(cntr).rho = rho;

        cntr = cntr + 1;
        
        cntrmu = cntrmu + 1;
    end
    cntrN = cntrN + 1;
    save(outfile3CDE, 'params', 'results');
end
params.mu = [];
params.N = [];
save(outfile3CDE, 'params', 'results');
