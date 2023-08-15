function [nhaplotypes, ngens, Zn, parentVector, mutationVector, bool_extinct] =  SimulateOneRealization_BranchingProcess(params)

p_geom = 1/(params.R0 + 1);     % get parameter p_geom for the geometric distribution

Zn = params.N0;         % this is the initial population size; scalar since all the same haplotype (wild-type); will develop into a matrix with rows = haplotypes, cols = generations
parentVector = NaN;     % column vector of length nhaplotypes that carries the parent identity of the haplotype
mutationVector = NaN;   % column vector of length nhaplotypes that carries the number of mutations that separates the haplotype from its parent

max_mutations = max(5, poissinv(0.99,params.mu));
mut_probabilities = poisspdf(0:max_mutations,params.mu); mut_probabilities(end) = mut_probabilities(end) + (1-sum(mut_probabilities));

bool_extinct = 0;
while 1
    [nTotHaplotypes, ngens] = size(Zn);
    Zn(1:nTotHaplotypes, ngens+1) = NaN;        % extend the Zn matrix out by one generation to place in new offspring from this generation
    
    tempZnMatrix = []; tempparentVector = []; tempmutationVector = [];
    for i = 1:nTotHaplotypes
        [offspringSameHaplotype, tempZnMatrix_thisHaplo, tempparentVector_thisHaplo, tempmutationVector_thisHaplo] = CreateOffspringOfOneHaplotype(Zn(i, ngens), i, ngens, p_geom, mut_probabilities);
        Zn(i, ngens+1) = offspringSameHaplotype;
        tempZnMatrix = [tempZnMatrix; tempZnMatrix_thisHaplo];
        tempparentVector = [tempparentVector; tempparentVector_thisHaplo];
        tempmutationVector = [tempmutationVector; tempmutationVector_thisHaplo];
    end
    Zn = [Zn; tempZnMatrix];
    parentVector = [parentVector; tempparentVector];
    mutationVector = [mutationVector; tempmutationVector];
    tot_Zn = sum(Zn(:,ngens+1));
    if tot_Zn > params.Nmax
        break; 
    end
    if tot_Zn == 0
        bool_extinct = 1;
        break; 
    end

end
[nhaplotypes, ngens] = size(Zn);
