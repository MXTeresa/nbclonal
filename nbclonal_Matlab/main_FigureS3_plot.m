function void = main_FigureS3_plot(void)

clear all; close all; clc;

params.R0 = 1.2;
params.mu = 0.2;
params.N = 2;

params.maxFinalSize = 80;
params.maxGeneratedMutantLineages = 200;
params.maxEstablishedMutantLineages = 200;
params.maxClonal = 6;

[n_clonals, rho] = calculateClonalPMF(params, 1)

infile = strcat('results_clonalpmf_N', int2str(params.N), '_R0', int2str(params.R0*10), '_mu', int2str(params.mu*10), '_lim', int2str(params.maxFinalSize), '_', int2str(params.maxGeneratedMutantLineages), '_', int2str(params.maxEstablishedMutantLineages));

load(infile)

params
results

x = 0:100;
y_overall = geopdf(x, params.pgeom);
y_wt = nbinpdf(x, params.negBin_r_wt, params.pgeom);
y_mut = nbinpdf(x, params.negBin_r_mut, params.pgeom);

figure(1);
subplot(2,3,1); bar(x, [y_overall; y_wt; y_mut]); 
xlabel('offspring number'); ylabel('probability');
legend('overall', 'wt', 'mutant'); axis([-0.9 8.9 0 1])

subplot(2,3,2); bar(results.finalSizes, results.prob_finalSize); 
xlabel('final size'); ylabel('probability'); 
yax = axis; axis([-0.9 50 0 yax(4)]);

subplot(2,3,3); bar(results.mutantLineagesGenerated, results.prob_mutantLineagesGenerated); 
xlabel('# mutant lineages generated'); ylabel('probability');
yax = axis; axis([-0.9 30 0 yax(4)]);

subplot(2,3,4); bar(results.mutantLineagesEstablished, results.prob_mutantLineagesEstablished);
xlabel('# mutant lineages established '); ylabel('probability');
yax = axis; axis([-0.9 10 0 yax(4)]);

P0_stacked = [results.S_Inf results.S_2plus];
PX_stacked = zeros(size(P0_stacked)); PX_stacked(1) = results.S_0;
P1plus_stacked = zeros(size(P0_stacked)); P1plus_stacked(1) = results.S_1;

subplot(2,3,5); bar([PX_stacked; P0_stacked; P1plus_stacked], 'stacked'); hold on; xlabel('outcome'); ylabel('probability');
xticks([1 2 3]); xticklabels({'P_X','P_0','P_{1+}'})
y = axis; axis([y(1) y(2) 0 1])
sum([PX_stacked P0_stacked P1plus_stacked])
%compare against analytical results:
plot([y(1) 1], [results.P_X_analytical results.P_X_analytical], 'k--');

subplot(2,3,6);
bar(0:params.maxClonal, results.rho); 
xlabel('# clonal variants (\rho)'); ylabel('probability');
y = axis; axis([y(1) y(2) 0 1])
