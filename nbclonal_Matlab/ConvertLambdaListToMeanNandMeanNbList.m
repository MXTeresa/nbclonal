function [meanNlist, meanNblist] = ConvertLambdaListToMeanNandMeanNbList(lambdalist, R0)

cntr = 1;
params.N0list = 1:50;
for lambda = lambdalist
    yvals = poisspdf(params.N0list, lambda);
    yvals_cond_unnorm = yvals.*(1-(1/R0).^params.N0list);
    yvals_cond_norm = yvals_cond_unnorm/sum(yvals_cond_unnorm);
    prob_N0list = yvals_cond_norm;
    meanNlist(cntr) = sum(params.N0list.*prob_N0list);
    for thisNb = 1:20
        prob_thisNb = 0;
        for thisN = params.N0list
            prob_thisNb = prob_thisNb + poisspdf(thisN, lambda)*binopdf(thisNb,thisN,1-(1/R0));
        end
        Nb_probs(thisNb) = prob_thisNb;
    end
    meanNblist(cntr) = sum((1:20).*Nb_probs/sum(Nb_probs));
    cntr = cntr + 1;
end
