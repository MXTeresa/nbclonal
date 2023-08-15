function void = main_simulateFigure3B_branchingProcess(void)

clear all; close all; clc;

params_emp.n_stoch_real = 100;           % number of stochastic simulations with successful invasion that we want, for each parameterization considered

params_emp.R0 = 1.6;                      % within-host R0 (= mean number of offspring of the geometric distribution)
params_emp.lambda = 2.1;                  % mean of Poisson distn for initial viral population size
params_emp.mu = 0.4;                     % per-genome, per infection cycle mutation rate (realized mutations = poisson random variable with mean mu)

params_emp.keepOnlySuccessfulInfections = 1;

params_emp.Nmax = 1e5;                  % simulation of the branching process stops once the viral population size has exceeded this threshold value 

n = 0;
    
outfile3A = strcat('figure3B_empirical_R0', int2str(params_emp.R0*10), '_mu', int2str(params_emp.mu*10), '_lambda', int2str(params_emp.lambda*10));
    
% simulate branching process model under a given parameterization until we get the overall number we want (n_stoch_real)
while n < params_emp.n_stoch_real
    n
    params_emp.N0 = poissrnd(params_emp.lambda);
    [nhaplotypes, ngens, Zn, parentVector, mutationVector, bool_extinct] =  SimulateOneRealization_BranchingProcess(params_emp);
    nhaplotypes_final = nhaplotypes;
    ngens_final = ngens;
    
    if bool_extinct 
        if ~params_emp.keepOnlySuccessfulInfections
            n = n + 1;
            results_emp.nclonal(n) = NaN;
            results_emp.extinct(n) = bool_extinct;
            results_emp
        end
    else
        [mutationList, maxMut] = GetMutationListHaplotypes(nhaplotypes_final, parentVector, mutationVector);
        [nclonal] = CalculateClonalMutations(Zn(:,ngens_final), nhaplotypes_final, mutationList, maxMut);
        
        n = n + 1;
        results_emp.nclonal(n) = nclonal;
        results_emp.extinct(n) = bool_extinct;
        save(outfile3A, 'params_emp', 'results_emp');
    end
    
    params.N0 = [];
    
end

save(outfile3A, 'params_emp', 'results_emp');
