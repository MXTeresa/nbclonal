function [finalSizes, prob_finalSize_N1, prob_finalSize_N] = GetFinalSizeWildtype(params, pi_wtExtinct1)

finalSizes = 1:params.maxFinalSize;

% calculate the final size distribution for a single particle:
p = NaN*ones(size(finalSizes));
p(1) = 1/((1+ (params.negBin_mean_wt/params.negBin_r_wt))^params.negBin_r_wt);
for y = 2:params.maxFinalSize
    j = 0:(y-2);
    p(y) = (prod(((j/params.negBin_r_wt) + y))/factorial(y))*((params.negBin_r_wt/(params.negBin_mean_wt+params.negBin_r_wt))^(params.negBin_r_wt*y))*((params.negBin_mean_wt*params.negBin_r_wt/(params.negBin_mean_wt+params.negBin_r_wt))^(y-1)); % equation 17 in Nishiura; equation (9) in Blumberg & Lloyd-Smith (2013) Inference
end

% leave off the dangling bit of probability mass.

% normalize this pmf, so that we are conditioning on wt extinction: 
prob_finalSize_N1 = p/pi_wtExtinct1; 

% numerically expand results to larger number of initial particles (N >= 1). This cannot be done analytically because pmf of
% is not a standard distribution.
prob_finalSize_N = ExpandFinalSizeFromOneN(params, prob_finalSize_N1, params.N); 
