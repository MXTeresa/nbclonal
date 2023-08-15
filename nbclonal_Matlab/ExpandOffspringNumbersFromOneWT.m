function p_N = ExpandOffspringNumbersFromOneWT(params, p, n_wt)

if n_wt == 1
    p_N = p;
    return;
else
    p_N_oneDown = ExpandOffspringNumbersFromOneWT(params, p, n_wt-1);
    p_N = zeros(size(0:params.maxGeneratedMutantLineages));
    for i = 0:params.maxGeneratedMutantLineages
        for j = 0:params.maxGeneratedMutantLineages
            if (i+j) <= params.maxGeneratedMutantLineages
                p_N(i+j+1) = p_N(i+j+1) + p(i+1)*p_N_oneDown(j+1);
            end
        end
    end
end