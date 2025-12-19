function T = build_table(Vgs_vec, Vds_vec, p)
% Returns struct with flattened arrays for all grid points
[VGS, VDS] = meshgrid(Vgs_vec, Vds_vec);
VGS = VGS(:); VDS = VDS(:);
n = numel(VGS);
Id  = zeros(n,1); gm = zeros(n,1); gds = zeros(n,1);
for i = 1:n
  [Id(i), gm(i), gds(i)] = ekv_gm_gds_fd(VGS(i), VDS(i), p);
end
T.VGS = VGS; T.VDS = VDS;
T.Id  = Id;  T.gm  = gm;  T.gds = gds;
end
