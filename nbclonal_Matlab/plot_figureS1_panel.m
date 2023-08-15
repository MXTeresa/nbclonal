function void = plot_figureS1_panel(params, results, results_emp, nrows, ncols, panelA, panelB, paneltitle)

results.logLMatrix_95 = results.logLMatrix;
locsMax = find(results.logLMatrix_95 == max(max(results.logLMatrix_95)));
maxLogL = results.logLMatrix_95(locsMax(1));
locsNaN = find(results.logLMatrix_95 < (maxLogL-2.995));
results.logLMatrix_95(locsNaN) = NaN;

subplot(nrows, ncols, panelB); s = pcolor(params.lambdalist, params.mulist, results.logLMatrix_95); s.LineStyle = 'none'; colorbar;
xlabel('\lambda'); ylabel('mutation rate \mu'); yax = axis; axis([0 4 0 4]); hold on;

maxLog = -Inf;
for i = 1:length(params.lambdalist)
    for j = 1:length(params.mulist)
        if results.logLMatrix(j,i) > maxLog
            bestLambdaLoc = i;
            bestMuLoc = j;
            maxLog = results.logLMatrix(j,i);
        end
    end
end
[params.lambdalist(bestLambdaLoc) params.mulist(bestMuLoc)]
plot([params.lambdalist(bestLambdaLoc) params.lambdalist(bestLambdaLoc)], [0 4], 'r--');
plot([0 10], [params.mulist(bestMuLoc) params.mulist(bestMuLoc)], 'r--');

clonal_list = 0:(max(results_emp.nclonal)+1);
for n_clonal = clonal_list
    n(n_clonal+1) = length(find(results_emp.nclonal == n_clonal));
end
freq_list_cond = (n/sum(n))';
n

% expectation of distribution under the MLE of lambda and mu:
params.mu = params.mulist(bestMuLoc); 
params.lambda = params.lambdalist(bestLambdaLoc);

params.maxFinalSize = 80;
params.maxGeneratedMutantLineages = 50;
params.maxEstablishedMutantLineages = 50;
params.maxClonal = max(clonal_list);

maxN = 8;
xvals = 0:maxN; yvals = poisspdf(xvals, params.lambda);
yvals_cond_unnorm = yvals.*(1-(1/params.R0).^xvals);
yvals_cond_norm_MLE = yvals_cond_unnorm/sum(yvals_cond_unnorm);
overallExpectedPMF = zeros(size(0:params.maxClonal));
for initPop = 0:maxN
    params.N = initPop;
    [n_clonals, rho] = calculateClonalPMF(params, 0)
    overallExpectedPMF = overallExpectedPMF + yvals_cond_norm_MLE(initPop+1)*rho;
end
[freq_list_cond overallExpectedPMF']
subplot(nrows, ncols, panelA); bar(clonal_list, [freq_list_cond overallExpectedPMF']);
xlabel('# clonal variants'); ylabel('proportion');
axis([-0.9 3.9 0 1.0]);
title(paneltitle)