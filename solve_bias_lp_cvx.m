function sol = solve_bias_lp_cvx(T, specs)
% struct with VDD, gm_min, eta_gmId, Av_min, Vds_min

n = numel(T.Id);
mask = (T.VDS >= specs.Vds_min) & (T.Id > 0) & (T.gm > 0) & (T.gds > 0);
Id  = T.Id(mask);  gm  = T.gm(mask);  gds = T.gds(mask);
VGS = T.VGS(mask); VDS = T.VDS(mask);
m = numel(Id);


cvx_solver sedumi
cvx_precision best


feas = (gm >= specs.gm_min) & ...
       (gm >= specs.eta_gmId .* Id) & ...
       (VDS >= specs.Vds_min);

if isfield(specs,'Av_min') && ~isempty(specs.Av_min)
    feas = feas & (gm >= specs.Av_min .* gds);
end

fprintf("Number of feasible grid points = %d\n", nnz(feas));
fprintf("Points meeting gm_min: %d\n", nnz(gm >= specs.gm_min));
fprintf("Points meeting eta:   %d\n", nnz(gm >= specs.eta_gmId .* Id));
fprintf("Points meeting Vds:   %d\n", nnz(VDS >= specs.Vds_min));
if isfield(specs,'Av_min') && ~isempty(specs.Av_min)
    fprintf("Points meeting Av:    %d\n", nnz(gm >= specs.Av_min .* gds));
end

% cvx_begin quiet
cvx_begin
  variable w(m)
  minimize( specs.VDD * (Id' * w) )
  subject to
    w >= 0;
    sum(w) == 1;

    (gm' * w) >= specs.gm_min;
    (gm' * w) >= specs.eta_gmId * (Id' * w);

    if isfield(specs,'Av_min') && ~isempty(specs.Av_min)
      (gm' * w) >= specs.Av_min * (gds' * w);
    end

    (VDS' * w) >= specs.Vds_min;
cvx_end

sol.cvx_status = cvx_status;
sol.optval = cvx_optval;

disp(cvx_status);
fprintf("optval = %.4e\n", cvx_optval);
fprintf("any NaN in w? %d\n", any(isnan(w)));
fprintf("sum(w)=%.6f, min(w)=%.3e\n", sum(w), min(w));
fprintf("gm(w)=%.3e, Id(w)=%.3e, gds(w)=%.3e, Vds(w)=%.3f\n", ...
    gm'*w, Id'*w, gds'*w, VDS'*w);

if any(isnan(w)) || ~strcmp(cvx_status,"Solved")
    error("Bad CVX solution. status=%s", cvx_status);
end




% [~,k] = max(w);

% Pick the MINIMUM current point among grid points that satisfy all constraints
feas = (gm >= specs.gm_min) & ...
       (gm >= specs.eta_gmId .* Id) & ...
       (VDS >= specs.Vds_min);

if isfield(specs,'Av_min') && ~isempty(specs.Av_min)
    feas = feas & (gm >= specs.Av_min .* gds);
end

if nnz(feas) == 0
    k = 1;  % fallback 
else
    Id_tmp = Id;
    Id_tmp(~feas) = Inf;   % exclude infeasible points
    [~,k] = min(Id_tmp);  % lowest current feasible point
end


sol.w = w;
sol.VGS_relaxed = VGS' * w;
sol.VDS_relaxed = VDS' * w;
sol.Id_relaxed  = Id'  * w;
sol.gm_relaxed  = gm'  * w;
sol.gds_relaxed = gds' * w;

sol.VGS_pick = VGS(k);
sol.VDS_pick = VDS(k);
end
