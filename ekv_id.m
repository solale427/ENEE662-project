function Id = ekv_id(Vgs, Vds, p)

% Simplified EKV expression 
argF = (p.kappa*(Vgs - p.VT0)) / (2*p.UT);
argR = (p.kappa*(Vgs - p.VT0) - Vds) / (2*p.UT);

F = log(1 + exp(argF));
R = log(1 + exp(argR));

Id0 = p.Is*(F.^2 - R.^2);
Id = Id0 .* (1 + p.lambda*max(Vds,0)); 
end
