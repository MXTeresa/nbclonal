function void = main_Figure2B_calculate(void)

clear all; close all; clc;

params_emp.R0list = 1.1:0.1:6;
params_emp.mulist = 0:0.001:0.1;
params_emp.N = 1;

params_emp.maxFinalSize = 120;
params_emp.maxGeneratedMutantLineages = 200;
params_emp.maxEstablishedMutantLineages = 200;
params_emp.maxClonal = 6;

outfile = 'figure2B';

mean_clonal_branchingProcess(1:length(params_emp.mulist),1:length(params_emp.R0list)) = NaN;

this_params_emp.N = params_emp.N;
this_params_emp.maxFinalSize = params_emp.maxFinalSize;
this_params_emp.maxGeneratedMutantLineages = params_emp.maxGeneratedMutantLineages;
this_params_emp.maxEstablishedMutantLineages = params_emp.maxEstablishedMutantLineages;
this_params_emp.maxClonal = params_emp.maxClonal;

cntr_mu = 1;
for this_mu = params_emp.mulist
    [this_mu length(params_emp.mulist)]
    cntr_R0 = 1;
    for this_R0 = params_emp.R0list
        this_params_emp.R0 = this_R0;
        this_params_emp.mu = this_mu;
        [n_clonals, rho] = calculateClonalPMF(this_params_emp, 0)

        mean_clonal_branchingProcess(cntr_mu,cntr_R0) = sum(n_clonals.*rho);
        cntr_R0 = cntr_R0 + 1;
    end
    cntr_mu = cntr_mu + 1;
    save(outfile, 'params_emp', 'mean_clonal_branchingProcess');
end

save(outfile, 'params_emp', 'mean_clonal_branchingProcess');
