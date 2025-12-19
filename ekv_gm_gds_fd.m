function [Id, gm, gds] = ekv_gm_gds_fd(Vgs, Vds, p)
dV = 1e-4;    % can adjust

Id  = ekv_id(Vgs, Vds, p);

Idp = ekv_id(Vgs + dV, Vds, p);
Idm = ekv_id(Vgs - dV, Vds, p);
gm  = (Idp - Idm) / (2*dV);

Idp = ekv_id(Vgs, Vds + dV, p);
Idm = ekv_id(Vgs, Vds - dV, p);
gds = (Idp - Idm) / (2*dV);
end
