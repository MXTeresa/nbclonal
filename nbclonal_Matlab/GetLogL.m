function logL = GetLogL(lambda, N0list, n_clonal, pmf_matrix, nclonal_data, R0)

yvals = poisspdf(N0list, lambda);
yvals_cond_unnorm = yvals.*(1-(1/R0).^N0list);
yvals_cond_norm = yvals_cond_unnorm/sum(yvals_cond_unnorm);
prob_N0list = yvals_cond_norm;

nTransmissionPairs = length(nclonal_data);
logL = 0;
for i = 1:nTransmissionPairs
    nClonal_emp = nclonal_data(i);
    loc = find(n_clonal == nClonal_emp);
    logL = logL + log(sum(pmf_matrix(:,loc).*(prob_N0list')));
end
