function void = main_Figure2_plot(void)

clear all; close all; clc;

params.R0list = 1.1:0.1:6;
params.mulist = 0:0.001:0.1;
params.N = 1;

bozic_delta = 1./params.R0list;
bozic_u = params.mulist;

cntr = 1;
for this_u = bozic_u
    bozic_mean_clonal(cntr,:) = bozic_delta*this_u./(1-bozic_delta);
    cntr = cntr + 1;
end
subplot(1,3,1); mesh(params.R0list, params.mulist, log10(bozic_mean_clonal)); xlabel('R_0'); ylabel('\mu'); zlabel('mean # of clonal variants');
xticks([1 3 5 7])
yticks([0 0.05 0.1])
zticks([-4 -2 0])
zticklabels({'10^{-4}','10^{-2}','10^0'})
xlim([1 7])
ylim([0 0.1])
zlim([-4 0])

clear all;

infile = 'figure2B';
load(infile);

subplot(1,3,2); mesh(params_emp.R0list, params_emp.mulist, log10(mean_clonal_branchingProcess)); xlabel('R_0'); ylabel('\mu'); zlabel('mean # of clonal variants');
xticks([1 3 5 7])
yticks([0 0.05 0.1])
zticks([-4 -2 0])
zticklabels({'10^{-4}','10^{-2}','10^0'})
xlim([1 7])
ylim([0 0.1])
zlim([-4 0])

clear all;

infile = 'figure2C_empirical_R012_mu2_N2';

load(infile);

n_emp_simulated = length(results_emp.nclonal)
clonal_list_emp = 0:(max(results_emp.nclonal)+1);
for n_clonal = clonal_list_emp
    n(n_clonal+1) = length(find(results_emp.nclonal == n_clonal));
end

empirical_clonalsims = n'
empirical_clonaltotalsims = sum(n)
pause

freq_list_cond_plot = (n/sum(n))';

params.R0 = params_emp.R0;
params.N = params_emp.N0;
params.mu = params_emp.mu;
params.maxFinalSize = 120;
params.maxGeneratedMutantLineages = 200;
params.maxEstablishedMutantLineages = 200;
params.maxClonal = max(6, max(clonal_list_emp));

[n_clonals, rho] = calculateClonalPMF(params, 1);

if max(clonal_list_emp) < 6
    freq_list_cond_plot_expanded = zeros(1,length(0:6));
    freq_list_cond_plot_expanded(1:(max(clonal_list_emp)+1)) = freq_list_cond_plot;
else
    freq_list_cond_plot_expanded = freq_list_cond_plot;
end

subplot(1,3,3); bar(n_clonals, [freq_list_cond_plot_expanded rho']);

xlabel('# clonal variants');
ylabel('frequency');
max_clonal_plot = 8.9; %max(clonal_list_emp) + 0.9;
axis([-0.9 max_clonal_plot 0 0.7])
