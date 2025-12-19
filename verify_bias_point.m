function out = verify_bias_point(Vgs, Vds, p, specs)
[Id, gm, gds] = ekv_gm_gds_fd(Vgs, Vds, p);

out.Vgs = Vgs; 
out.Vds = Vds;
out.Id = Id; 
out.gm = gm; 
out.gds = gds;
out.P  = specs.VDD * Id;
out.gmId = gm / Id;
out.Av = gm / gds;

out.pass_gm   = (gm >= specs.gm_min);
out.pass_eta  = (gm >= specs.eta_gmId * Id);
out.pass_Vds  = (Vds >= specs.Vds_min);
if isfield(specs,'Av_min') && ~isempty(specs.Av_min)
    out.pass_Av = (out.Av >= specs.Av_min);
else
   out.pass_Av = true;
end
  out.pass_all = out.pass_gm && out.pass_eta && out.pass_Vds && out.pass_Av;
end
