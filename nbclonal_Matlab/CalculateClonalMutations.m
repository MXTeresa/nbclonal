function nclonal = CalculateClonalMutations(finalPopSizes, nhaplotypes_final, mutationList, maxMut)

totFinalPopSize = sum(finalPopSizes);

mutFinalPopSizes = zeros(maxMut, 1);

for j = 1:nhaplotypes_final
    this_list = mutationList(j).muts;
    for k = this_list
        mutFinalPopSizes(k,1) = mutFinalPopSizes(k,1) + finalPopSizes(j);
    end
end

locs_clonal = find(mutFinalPopSizes == totFinalPopSize);
nclonal = length(locs_clonal); % clonal mutations are ones that the whole population carries
