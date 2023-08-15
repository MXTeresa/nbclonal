function pi_Extinct = GetProbabilityOfExtinction(negBin_mean, negBin_r)

if negBin_mean <= 1
    pi_Extinct = 1;
    return;
end

if isinf(negBin_r)
    %pgf = exp(-R0*(1-s)) % constant R0, p. 12 JLS supp
    getProbExt = @(x)(abs(x-exp(-negBin_mean*(1-x))));
else
    %pgf = (1+(R0/k)*(1-s))^(-k); %   JLS Methods section
    getProbExt = @(x)(abs(x-(1+(negBin_mean/negBin_r)*(1-x))^(-negBin_r))); 
end

x = fminbnd(getProbExt,0,1);
pi_Extinct = x;
    
% this is the transcendental function:
%pi_Extinct = 1/(1+((negBin_mean*(1-pi_Extinct))/negBin_r)^negBin_r;   % equation (4) Nishiura; also JLS paper
    