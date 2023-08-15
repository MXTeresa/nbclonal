function void = main_Figure3_plot(void)

clear all; close all; clc;

this_lambda = 2.1; this_mu = 0.4; this_R0 = 1.6;
xvals = 0:10; yvals = poisspdf(xvals, this_lambda);
yvals_cond_unnorm = yvals.*(1-(1/this_R0).^xvals);
yvals_cond_norm = yvals_cond_unnorm/sum(yvals_cond_unnorm);

subplot(3,2,1); bar(xvals, [yvals' yvals_cond_norm']);
xlabel('initial viral population size N'); ylabel('probability');
axis([-0.9 8.9 0 0.4])

infile = 'figure3B_empirical_R016_mu4_lambda21'; load(infile);
load(infile);
clonal_list = 0:(max(results_emp.nclonal)+1);
for n_clonal = clonal_list
    n(n_clonal+1) = length(find(results_emp.nclonal == n_clonal));
end

freq_list_cond = (n/params_emp.n_stoch_real)';
n_clonal'
n'
sum(n)

% expectation of distribution under the MLE of lambda and mu:

lambdaMLE = 2.3;
muMLE = 0.39;

params.R0 = 1.6;
params.mu = muMLE;
params.lambda = lambdaMLE;

params.maxFinalSize = 120;
params.maxGeneratedMutantLineages = 200;
params.maxEstablishedMutantLineages = 200;
params.maxClonal = max(clonal_list);

xvals = 0:10; yvals = poisspdf(xvals, params.lambda);
yvals_cond_unnorm = yvals.*(1-(1/params.R0).^xvals);
yvals_cond_norm_MLE = yvals_cond_unnorm/sum(yvals_cond_unnorm);

overallExpectedPMF = zeros(size(0:8));
for initPop = 0:8
    params.N = initPop;
    [n_clonals, rho] = calculateClonalPMF(params, 0)
    overallExpectedPMF = overallExpectedPMF + yvals_cond_norm_MLE(initPop+1)*rho;
end

subplot(3,2,2); bar(clonal_list, [freq_list_cond overallExpectedPMF']);

xlabel('# clonal variants'); ylabel('proportion');
axis([-0.9 8.9 0 0.85])

clear all;

load('figure3CDE_data')

cntr = 1;
cntrN = 1;
for N0 = params.N0list
    params.N = N0;
    cntrmu = 1;
    for mu = params.mulist
        [N0 mu]
        params.mu = mu;
        results.muNvals(cntr,1:2) = [N0 mu];

        results.prob0(cntrmu, cntrN) = results.rhovals(cntr).rho(1);
        results.prob1(cntrmu, cntrN) = results.rhovals(cntr).rho(2);
        results.prob2(cntrmu, cntrN) = results.rhovals(cntr).rho(3);
        
        cntr = cntr + 1;
        cntrmu = cntrmu + 1;
    end
    cntrN = cntrN + 1;
end

loc = max(find(params.mulist <= 1.5));

subplot(3,3,4); s = pcolor(params.N0list, params.mulist(1:loc), results.prob0(1:loc,:)); s.LineStyle = 'none'; colorbar; xlabel('N'); ylabel('\mu'); 
xticks([0.5:8.5]); headers = ["0", "1","2", "3", "4", "5","6","7","8"]; xticklabels(headers);
%title('prob. of observing 0 clonal variants');
subplot(3,3,5); s = pcolor(params.N0list, params.mulist(1:loc), results.prob1(1:loc,:)); s.LineStyle = 'none'; colorbar; xlabel('N'); ylabel('\mu'); 
xticks([0.5:8.5]); headers = ["0", "1","2", "3", "4", "5","6","7","8"]; xticklabels(headers);
%title('prob. of observing 1 clonal variant');
subplot(3,3,6); s = pcolor(params.N0list, params.mulist(1:loc), results.prob2(1:loc,:)); s.LineStyle = 'none'; colorbar; xlabel('N'); ylabel('\mu'); 
%title('prob. of observing 2 clonal variants');
xticks([0.5:8.5]); headers = ["0", "1","2", "3", "4", "5","6","7","8"]; xticklabels(headers);


clear all; clc;

infile = 'figure3FGH_data';
load(infile);

results.logLMatrix_95 = results.logLMatrix;
locsMax = find(results.logLMatrix_95 == max(max(results.logLMatrix_95)));
maxLogL = results.logLMatrix_95(locsMax(1));
locsNaN = find(results.logLMatrix_95 < (maxLogL-2.995));
results.logLMatrix_95(locsNaN) = NaN;

subplot(3,3,7); s = pcolor(params.lambdalist, params.mulist, results.logLMatrix_95); s.LineStyle = 'none'; colorbar;
xlabel('\lambda'); ylabel('mutation rate \mu'); yax = axis; axis([0 10 0 1.5]); hold on;

params.lambda = 2.1;
plot([params.lambda params.lambda], [0 2], 'k--');

plot([0 10], [0.4 0.4], 'k--');

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
plot([params.lambdalist(bestLambdaLoc) params.lambdalist(bestLambdaLoc)], [0 1.5], 'r--');
plot([0 10], [params.mulist(bestMuLoc) params.mulist(bestMuLoc)], 'r--');

subplot(3,3,8); s = pcolor(params.meanNlist, params.mulist, results.logLMatrix_95); s.LineStyle = 'none'; colorbar;
xlabel('mean N'); ylabel('mutation rate \mu'); yax = axis; axis([0 10 0 1.5])

subplot(3,3,9); s = pcolor(params.meanNblist, params.mulist, results.logLMatrix_95); s.LineStyle = 'none'; colorbar;
xlabel('mean N_b'); ylabel('mutation rate \mu'); yax = axis; axis([0 10 0 1.5])

