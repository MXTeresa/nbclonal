function prob_finalSize = ExpandFinalSizeFromOneN(params, prob_finalSize_N1, N)

p = prob_finalSize_N1;

if N == 1
    prob_finalSize = p;
    return;
else
    p_N_oneDown = ExpandFinalSizeFromOneN(params, p, N-1);
    p_N = zeros(size(1:params.maxFinalSize));
    for i = 1:params.maxFinalSize
        for j = 1:params.maxFinalSize
            if (i+j) <= params.maxFinalSize
                p_N(i+j) = p_N(i+j) + p(i)*p_N_oneDown(j);
            end
        end
    end
end

% leave off the dangling bit of probability mass
prob_finalSize = p_N;
