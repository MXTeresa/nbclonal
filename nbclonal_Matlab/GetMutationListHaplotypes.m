function [mutationList, maxMut] = GetMutationListHaplotypes(nhaplotypes_final, parentVector, mutationVector)

% returns a list of mutations that each haplotype carries

mutationList(1).muts = []; % the original haplotype carries no mutations
maxMut = 0;
for i = 2:nhaplotypes_final
    nMuts = mutationVector(i);
    extraMuts = maxMut + [1:nMuts];
    mutationList(i).muts = [mutationList(parentVector(i)).muts extraMuts];
    maxMut = max(extraMuts);
end
