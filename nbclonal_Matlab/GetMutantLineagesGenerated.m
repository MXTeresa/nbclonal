function [mutantLineagesGenerated, prob_mutantLineagesGenerated] = GetMutantLineagesGenerated(params, finalSizes, prob_finalSize)

mutantLineagesGenerated = 0:params.maxGeneratedMutantLineages;

if params.mu == 0
    prob_mutantLineagesGenerated = zeros(size(mutantLineagesGenerated));
    prob_mutantLineagesGenerated(1) = 1;
    return;
end

% calculate the pmf for the number of mutant lineages generated, given final size distn pmf provided.

prob_mutantLineagesGenerated = zeros(size(mutantLineagesGenerated));
prob_exactly_this_many_mut_offspring = nbinpdf(mutantLineagesGenerated, params.negBin_r_mut, params.pgeom); % this is from a single wt particle

for n_wt = finalSizes
    % this sometimes runs into numerical issues:
    %prob_mutantLineagesGenerated = prob_mutantLineagesGenerated + prob_finalSize(n_wt)*nbinpdf(mutantLineagesGenerated, n_wt*params.negBin_r_mut, params.pgeom);
    % so use this code instead:
    prob_exactly_this_many_mut_offspring_N = ExpandOffspringNumbersFromOneWT(params, prob_exactly_this_many_mut_offspring, n_wt);
    prob_mutantLineagesGenerated = prob_mutantLineagesGenerated + prob_finalSize(n_wt)*prob_exactly_this_many_mut_offspring_N;
end
