clear; clc; close all;

p = default_params();

% ---- Base specs ----
specs.VDD = 1.2;
specs.gm_min = 1e-3;        % 1 mS
specs.Av_min = 20;          % keep this or set [] to disable
specs.Vds_min = 0.2;

% We'll sweep this:
% specs.eta_gmId = ???;

Vgs_vec = linspace(0.4,0.9,101);
Vds_vec = linspace(0.1,1.2,111);

% ---- Build table and auto-scale Is once so gm_min is achievable somewhere ----
T = build_table(Vgs_vec, Vds_vec, p);
gm_max = max(T.gm);
fprintf("gm_max = %.3e S\n", gm_max);

if gm_max < specs.gm_min
    scale = (specs.gm_min / gm_max) * 1.2;   % 20% margin
    fprintf("Scaling Is by %.2f to meet gm_min...\n", scale);
    p.Is = p.Is * scale;
    T = build_table(Vgs_vec, Vds_vec, p);
end

fprintf("gm_max after scaling = %.3e S\n", max(T.gm));

% ==========================================================
% ETA SWEEP: tradeoff curve + feasibility boundary
% ==========================================================
etas = 2:1:20;              % sweep range (adjust later if you want)
Pvals = NaN(size(etas));    % power (W)
Idvals = NaN(size(etas));   % current (A)
VGSvals = NaN(size(etas));  % chosen VGS (V)
VDSvals = NaN(size(etas));  % chosen VDS (V)
feasible = false(size(etas));

for t = 1:numel(etas)
    specs.eta_gmId = etas(t);

    try
        sol = solve_bias_lp_cvx(T, specs);

        if strcmp(sol.cvx_status,"Solved")
            feasible(t) = true;

            ver = verify_bias_point(sol.VGS_pick, sol.VDS_pick, p, specs);

            Pvals(t)   = ver.P;
            Idvals(t)  = ver.Id;
            VGSvals(t) = ver.Vgs;
            VDSvals(t) = ver.Vds;
        else
            feasible(t) = false;
            Pvals(t)   = NaN;
            Idvals(t)  = NaN;
            VGSvals(t) = NaN;
            VDSvals(t) = NaN;
        end

    catch
        feasible(t) = false;
        Pvals(t)   = NaN;
        Idvals(t)  = NaN;
        VGSvals(t) = NaN;
        VDSvals(t) = NaN;
    end
end


% ---- Print feasibility boundary ----
idx = find(feasible);
if isempty(idx)
    disp("No feasible eta values found in this sweep range.");
else
    fprintf("Feasible eta range: [%g, %g]\n", etas(idx(1)), etas(idx(end)));
end

% ---- Plot 1: Power vs eta ----
figure;
P_mW = 1e3 * Pvals;                 % W -> mW
plot(etas, P_mW, '-o'); grid on;
xlabel('\eta = (g_m/I_D)_{min}  [1/V]');
ylabel('Power P = V_{DD} I_D  [mW]');
title('Tradeoff: Minimum Efficiency Constraint vs Required Power');

mask = isfinite(P_mW);
if any(mask)
    ylim([0, 1.2*max(P_mW(mask))]);
end


% ---- Plot 2: Current vs eta ----
figure;
Id_mA = 1e3 * Idvals;               % A -> mA
plot(etas, Id_mA, '-o'); grid on;
xlabel('\eta = (g_m/I_D)_{min}  [1/V]');
ylabel('Drain current I_D  [mA]');
title('Tradeoff: Minimum Efficiency Constraint vs Required Current');

mask = isfinite(Id_mA);
if any(mask)
    ylim([0, 1.2*max(Id_mA(mask))]);
end

% ---- (Optional) Plot 3: Chosen bias voltages vs eta ----
figure;
plot(etas, VGSvals, '-o'); hold on;
plot(etas, VDSvals, '-o'); grid on;
xlabel('\eta = (g_m/I_D)_{min}  [1/V]');
ylabel('Voltage [V]');
legend('V_{GS}^*','V_{DS}^*','Location','best');
title('Optimal Bias Voltages vs Efficiency Constraint');

mask = isfinite([VGSvals(:); VDSvals(:)]);
if any(mask)
    lo = min([VGSvals(:); VDSvals(:)],[],'omitnan');
    hi = max([VGSvals(:); VDSvals(:)],[],'omitnan');
    ylim([lo-0.05, hi+0.05]);
end


% ==========================================================
% OPTIONAL: pick one eta and show a final single solution
% ==========================================================
% Choose a feasible eta near the boundary (most interesting):
if ~isempty(idx)
    specs.eta_gmId = etas(idx(end));   % tightest feasible eta
    sol = solve_bias_lp_cvx(T, specs);
    ver = verify_bias_point(sol.VGS_pick, sol.VDS_pick, p, specs);

    disp("Final chosen eta solution (tight feasible):");
    disp(ver);
end


figure;
stem(etas, feasible, 'filled'); grid on;
xlabel('\eta = (g_m/I_D)_{min} [1/V]');
ylabel('Feasible? (1=yes, 0=no)');
title('Feasibility vs Efficiency Constraint');
ylim([-0.1 1.1]);
