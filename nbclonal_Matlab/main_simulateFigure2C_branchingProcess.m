function void = main_simulateFigure2C_branchingProcess(void)

clear all; close all; clc;

rng('shuffle')

% Figure 2C:
params_emp.R0 = 1.2;            % within-host R0 (= mean number of offspring)
params_emp.N0 = 2;              % initial population size
params_emp.mu = 0.2;            % mean number of mutations that occur at replication (realized mutations = poisson random variable with mean mu)

params_emp.n_stoch_real = 4000;     % the number of stochastic realizations used to generate the pmf of clonal and subclonal mutations
params_emp.Nmax = 1e5;              % simulation of the branching process stops once the population size has exceeded this threshold value 

params_emp.keepOnlySuccessfulInfections = 1;

outfile2C = strcat('figure2C_empirical_R0', int2str(round(params_emp.R0*10)), '_mu', int2str(params_emp.mu*10), '_N', int2str(params_emp.N0));

n = 0;

% branching process simulations:
while n < params_emp.n_stoch_real
    n
    
    [nhaplotypes, ngens, Zn, parentVector, mutationVector, bool_extinct] =  SimulateOneRealization_BranchingProcess(params_emp);
    nhaplotypes_final = nhaplotypes;
    ngens_final = ngens;
    
    if bool_extinct 
        if ~params_emp.keepOnlySuccessfulInfections
            n = n + 1;
            results_emp.nclonal(n) = NaN;
            results_emp.extinct(n) = bool_extinct;
        end
    else
        [mutationList, maxMut] = GetMutationListHaplotypes(nhaplotypes_final, parentVector, mutationVector);
        [nclonal] = CalculateClonalMutations(Zn(:,ngens_final), nhaplotypes_final, mutationList, maxMut);
        
        n = n + 1;
        results_emp.nclonal(n) = nclonal;
        results_emp.extinct(n) = bool_extinct;
        save(outfile2C, 'params_emp', 'results_emp');
    end
    
end

save(outfile2C, 'params_emp', 'results_emp');
