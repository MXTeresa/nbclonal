function [offspringSameHaplotype, tempZnMatrix, tempparentVector, tempmutationVector] = CreateOffspringOfOneHaplotype(n_thisHaplotype, hap_num, ngens, p_geom, mut_probabilities)

Znplus1_thisHaplotype = geornd(p_geom, [n_thisHaplotype 1]);
tot_Znplus1 = sum(Znplus1_thisHaplotype);       % gets total offspring from that haplotype

tempZnMatrix = []; tempparentVector = []; tempmutationVector = [];
    
if tot_Znplus1 == 0
    offspringSameHaplotype = 0;
    return;
end

haplotypeOffspring = mnrnd(tot_Znplus1, mut_probabilities, 1);

offspringSameHaplotype = haplotypeOffspring(1);
for j = 2:length(mut_probabilities)
    if haplotypeOffspring(j) > 0
        n_mutations = j-1;
        insertMatrix = [zeros(haplotypeOffspring(j), length(1:ngens)) ones(haplotypeOffspring(j),1)];
        tempZnMatrix = [tempZnMatrix; insertMatrix];
        tempparentVector = [tempparentVector; ones(haplotypeOffspring(j),1)*hap_num];
        tempmutationVector = [tempmutationVector; ones(haplotypeOffspring(j),1)*n_mutations]; % number of mutations from parent
    end
end
