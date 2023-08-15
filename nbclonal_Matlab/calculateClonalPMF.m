function [n_clonals, rho] = calculateClonalPMF(params, save_bool)

if params.N < 1
    display('N has to be an integer that is greater than or equal to 1');
    display('returning 0 values for rho')
    n_clonals = 0:params.maxClonal;
    rho = zeros(size(n_clonals));
    return;
end

outfile = strcat('results_clonalpmf_N', int2str(params.N), '_R0', int2str(params.R0*10), '_mu', int2str(params.mu*10), '_lim', int2str(params.maxFinalSize), '_', int2str(params.maxGeneratedMutantLineages), '_', int2str(params.maxEstablishedMutantLineages));

% parameterize the overall geometric offspring distribution, the wt offspring distribution, 
% and the mutant offspring distribution from wt particles:
params.pgeom = 1/(params.R0 + 1);
params.negBin_r_wt = exp(-params.mu);
params.negBin_r_mut = 1-exp(-params.mu);

params.negBin_mean_wt = params.negBin_r_wt*params.R0;
params.negBin_mean_mut = params.negBin_r_mut*params.R0;

% calculate the probability of extinction of the overall viral population (this is S_0, equiv. to P_X):
if params.R0 <= 1
    results.P_X_analytical = 1;
else
    results.P_X_analytical = (1/params.R0)^params.N;
end

% calculate the probability of the wild-type viral population establishing (s_inf, contributes to p_0):
% to do this, first calculate the probability that a single wt particle's wt lineage goes extinct:
pi_wtExtinct1 = GetProbabilityOfExtinction(params.negBin_mean_wt, params.negBin_r_wt);
% then, calculate the probability that the overall wt population, of initial size N, establishes:
results.S_Inf = 1-pi_wtExtinct1^params.N; 

% calculate the final size distribution of minor outbreaks for wild-type
% viral particles (conditional on wild-type extinction, such that the pmf sums to 1):
[finalSizes, prob_finalSize_N1, prob_finalSize_N] = GetFinalSizeWildtype(params, pi_wtExtinct1);
results.finalSizes = finalSizes;
results.prob_finalSize_N1 = prob_finalSize_N1;
results.prob_finalSize = prob_finalSize_N;

% calculate pmf of the number of mutant lineages generated, given final size distribution pmf from above:
[mutantLineagesGenerated, prob_mutantLineagesGenerated] = GetMutantLineagesGenerated(params, results.finalSizes, results.prob_finalSize);
results.mutantLineagesGenerated = mutantLineagesGenerated;
results.prob_mutantLineagesGenerated = prob_mutantLineagesGenerated;

% calculate pmf for the number of mutant lineages that *establish*:
[mutantLineagesEstablished, prob_mutantLineagesEstablished] = GetMutantLineageEstablished(params, results.prob_mutantLineagesGenerated);
results.mutantLineagesEstablished = mutantLineagesEstablished;
results.prob_mutantLineagesEstablished = prob_mutantLineagesEstablished;

% use these probabilities to get S_0, S_1, S_2, etc.:
results.S_0 = (1-results.S_Inf)*results.prob_mutantLineagesEstablished(find(results.mutantLineagesEstablished == 0));
results.S_1 = (1-results.S_Inf)*results.prob_mutantLineagesEstablished(find(results.mutantLineagesEstablished == 1));
results.S_2plus = (1-results.S_Inf)*results.prob_mutantLineagesEstablished(find(results.mutantLineagesEstablished >= 2));

%display('S_{inf}, S_0, S_1, S_{2+} probabilities should add up to 1. In practice, will be slightly lower because we have to set a cutoff n for S_n.')
%sum([results.S_Inf results.S_0 results.S_1 results.S_2plus])

results.P_X = results.S_0; 
results.P_0_coarse = results.S_Inf + sum(results.S_2plus);
results.P_1plus = results.S_1;
results.P_0 = 1-results.S_0-results.S_1;

if results.P_X_analytical == 1
    results.rho0 = NaN;
    results.rho1plus = NaN;
    results.n_clonals = 0:params.maxClonal;
    results.P_resolved_clonals = NaN*ones(size(results.n_clonals));
else
    [n_clonals, P_resolved] = ResolveClonals(params, results.P_X, results.P_0, results.P_1plus);
    results.n_clonals = n_clonals;
    results.P_resolved_clonals = P_resolved;

    results.rho0 = results.P_0/(1-results.P_X_analytical);
    results.rho1plus = results.P_1plus/(1-results.P_X_analytical);
end

results.rho = results.P_resolved_clonals/(1-results.P_X_analytical);
nonAdjusted_rho = results.rho;
results.rho = results.rho/sum(results.rho); % should just readjust if some distributions come up short; user should check that they don't
adjusted_rho = results.rho;

n_clonals = results.n_clonals;
rho = results.rho;

if save_bool
    save(outfile, 'params', 'results');
end
