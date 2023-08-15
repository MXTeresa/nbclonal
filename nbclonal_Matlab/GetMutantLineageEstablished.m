function [mutantLineagesEstablished, prob_mutantLineagesEstablished] = GetMutantLineageEstablished(params, prob_mutantLineagesGenerated)

% Use a binomial distribution with success probability of 1-1/R0 (given the mutant's geometric offspring distribution, and M trials, where M is the number of mutant lineage offspring

mutantLineagesEstablished = 0:params.maxEstablishedMutantLineages;

prob_mutantLineagesEstablished = [];
for i = mutantLineagesEstablished
    binoVals = [];
    for M = 0:params.maxGeneratedMutantLineages
        if params.R0 <= 1
            binoVals = [binoVals binopdf(i,M,0)];
        else
            binoVals = [binoVals binopdf(i,M,1-1/params.R0)];
        end
    end
    prob_mutantLineagesEstablished = [prob_mutantLineagesEstablished sum(prob_mutantLineagesGenerated.*binoVals)];
end
prob_mutantLineagesEstablished(end) = prob_mutantLineagesEstablished(end) + (1-sum(prob_mutantLineagesEstablished));
