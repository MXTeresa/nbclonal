function [n_clonals, P_resolved] = ResolveClonals(params, P_X, P_0, P_1plus)

n_clonals = 0:params.maxClonal;

if P_1plus == 0
    P_resolved = zeros(size(n_clonals));
    P_resolved(1) = P_0;
    return;
end

%P_resolved(1) = P_0;  % 0 clonals

% figure out P_0, P_1, P_2, etc. for n = 1 case first:
paramsN1 = params; paramsN1.N = 1;
%resultsN1.P_X_analytical = (1/params.R0)^paramsN1.N;
pi_wtExtinct1 = GetProbabilityOfExtinction(paramsN1.negBin_mean_wt, paramsN1.negBin_r_wt);
% then, calculate the probability that the overall wt population, of initial size N, establishes:
resultsN1.S_Inf = 1-pi_wtExtinct1^paramsN1.N; 
[finalSizes, prob_finalSize_N1, prob_finalSize_N] = GetFinalSizeWildtype(paramsN1, pi_wtExtinct1);
resultsN1.finalSizes = finalSizes;
resultsN1.prob_finalSize_N1 = prob_finalSize_N1;
resultsN1.prob_finalSize = prob_finalSize_N;
[mutantLineagesGenerated, prob_mutantLineagesGenerated] = GetMutantLineagesGenerated(paramsN1, resultsN1.finalSizes, resultsN1.prob_finalSize);
resultsN1.mutantLineagesGenerated = mutantLineagesGenerated;
resultsN1.prob_mutantLineagesGenerated = prob_mutantLineagesGenerated;
[mutantLineagesEstablished, prob_mutantLineagesEstablished] = GetMutantLineageEstablished(paramsN1, resultsN1.prob_mutantLineagesGenerated);
resultsN1.mutantLineagesEstablished = mutantLineagesEstablished;
resultsN1.prob_mutantLineagesEstablished = prob_mutantLineagesEstablished;

resultsN1.S_0 = (1-resultsN1.S_Inf)*resultsN1.prob_mutantLineagesEstablished(find(resultsN1.mutantLineagesEstablished == 0));
resultsN1.S_1 = (1-resultsN1.S_Inf)*resultsN1.prob_mutantLineagesEstablished(find(resultsN1.mutantLineagesEstablished == 1));
resultsN1.P_0 = 1-resultsN1.S_0-resultsN1.S_1;
resultsN1.P_X = 1/params.R0;
resultsN1.P_1plus = resultsN1.S_1;

N1_Plist = NaN*ones(size(n_clonals));
k = 0;
N1_Plist(k+1) = resultsN1.P_0;

for k = 1:params.maxClonal
    sumioverk = 0;
    for i = 1:k
        index_val = k-i;
        sumioverk = sumioverk + (poisspdf(i,paramsN1.mu)/(1-poisspdf(0,paramsN1.mu)))*(N1_Plist(index_val+1)/(1-resultsN1.P_X));
    end
    N1_Plist(k+1) = resultsN1.P_1plus*sumioverk;
end
%N1_Plist
%[resultsN1.P_1plus sum(N1_Plist(2:end))]

N_Plist = NaN*ones(size(n_clonals));
k = 0;
N_Plist(k+1) = P_0;

for k = 1:params.maxClonal
    sumioverk = 0;
    for i = 1:k
        index_val = k-i;
        sumioverk = sumioverk + (poisspdf(i,params.mu)/(1-poisspdf(0,params.mu)))*(N1_Plist(index_val+1)/(1-resultsN1.P_X));
    end
    N_Plist(k+1) = P_1plus*sumioverk;
end
%N_Plist
%[P_1plus sum(N_Plist(2:end))]

%length(n_clonals)
%length(N_Plist)

P_resolved = N_Plist;

%P_resolved(2) = P_1plus*(poisspdf(1,params.mu)/(1-poisspdf(0,params.mu)))*(resultsN1.P_0/(1-resultsN1.P_X)); % exactly one clonal
%P_resolved_approx = P_1plus*(poisspdf(1,params.mu)/(1-poisspdf(0,params.mu)));

%[P_X, P_0, P_1plus P_resolved_approx P_resolved(2)]

%pause
