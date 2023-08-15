function void = main_Figure3FGH_calculate(void)

clear all; close all; clc;

infile = 'figure3B_empirical_R016_mu4_lambda21'; load(infile);

load(infile);
emp_clonal_list = 0:max(results_emp.nclonal);
for n_clonal = emp_clonal_list
    emp_n(n_clonal+1) = length(find(results_emp.nclonal == n_clonal));
end

params.R0 = 1.6;
params.mulist = 0.01:0.01:2;
params.lambdalist = 0.1:0.1:8;

% since can't go up to Inf, use these max values:
params.maxFinalSize = 120;
params.maxGeneratedMutantLineages = 200;
params.maxEstablishedMutantLineages = 200;
params.maxClonal = 10;

outfile3FGH = 'figure3FGH_data';

cntrmu = 1;
for mu = params.mulist
    mu
    params.mu = mu;
    N0_matrix = [];
    cntrN = 1; params.N0list = 1:20;
    for N0 = params.N0list
        params.N = N0;
        [n_clonals, rho] = calculateClonalPMF(params, 0)
        pmf_matrix(cntrN,:) = rho;
        cntrN = cntrN + 1;
    end
    cntrlambda = 1;
    for lambda = params.lambdalist
        logL_val = GetLogL(lambda, params.N0list, n_clonals, pmf_matrix, results_emp.nclonal, params.R0);
        logL_val
        results.logLMatrix(cntrmu, cntrlambda) = logL_val;
        cntrlambda = cntrlambda + 1;
    end
    cntrmu = cntrmu + 1;
end

save(outfile3FGH, 'params', 'results');

[params.meanNlist, params.meanNblist] = ConvertLambdaListToMeanNandMeanNbList(params.lambdalist, params.R0);

save(outfile3FGH, 'params', 'results');

