clear; clc; close all;
p = default_params();

% Specs
specs.VDD      = 1.2;
specs.gm_min   = 1e-3;   % 1 mS
specs.eta_gmId = 4;      % 1/V
specs.Av_min   = 20;     % intrinsic gain min
specs.Vds_min  = 0.2;    % V

Vgs_vec = linspace(0.4, 0.9, 101);
Vds_vec = linspace(0.1, 1.2, 111);

T = build_table(Vgs_vec, Vds_vec, p);

gm_max = max(T.gm);
fprintf("gm_max = %.3e S\n", gm_max);

if gm_max < specs.gm_min
    scale = (specs.gm_min / gm_max) * 1.2;
    fprintf("Scaling Is by %.2f to meet gm_min\n", scale);
    p.Is = p.Is * scale;
    T = build_table(Vgs_vec, Vds_vec, p);
end
fprintf("gm_max after scaling = %.3e S\n", max(T.gm));

%  Solve optimization
sol = solve_bias_lp_cvx(T, specs);
ver = verify_bias_point(sol.VGS_pick, sol.VDS_pick, p, specs);

disp(sol);
disp(ver);

VGSopt = sol.VGS_pick;
VDSopt = sol.VDS_pick;

gmId = T.gm ./ T.Id;
Av   = T.gm ./ T.gds;

[VGSg, VDSg] = meshgrid(Vgs_vec, Vds_vec);
IDg   = reshape(T.Id,  numel(Vds_vec), numel(Vgs_vec));
GMg   = reshape(T.gm,  numel(Vds_vec), numel(Vgs_vec));
GDSg  = reshape(T.gds, numel(Vds_vec), numel(Vgs_vec));
GMIDg = reshape(gmId,  numel(Vds_vec), numel(Vgs_vec));
AVg   = reshape(Av,    numel(Vds_vec), numel(Vgs_vec));


% EKV plots
VDS_slices = [0.2 0.6 1.2];
VGS_slices = [0.60 0.75 0.85];

% id vs vgs
figure; hold on; grid on;
for k = 1:numel(VDS_slices)
    [~, j] = min(abs(Vds_vec - VDS_slices(k)));
    plot(Vgs_vec, IDg(j,:), 'LineWidth', 1.5);
end
xlabel('V_{GS} [V]');
ylabel('I_D [mA]');
title('EKV sanity: I_D vs V_{GS} (fixed V_{DS})');
legend(compose('V_{DS}=%.1f V', VDS_slices));

% id vds
figure; hold on; grid on;
for k = 1:numel(VGS_slices)
    [~, i] = min(abs(Vgs_vec - VGS_slices(k)));
    plot(Vds_vec, IDg(:,i), 'LineWidth', 1.5);
end
xlabel('V_{DS} [V]');
ylabel('I_D [mA]');
title('EKV sanity: I_D vs V_{DS} (fixed V_{GS})');
legend(compose('V_{GS}=%.2f V', VGS_slices), 'Location','best');

% heatmaps
plotHeat(Vgs_vec, Vds_vec, GMg,   'g_m [S]',         'Heatmap: g_m(V_{GS},V_{DS})', VGSopt, VDSopt);
plotHeat(Vgs_vec, Vds_vec, GDSg,  'g_{ds} [S]',      'Heatmap: g_{ds}(V_{GS},V_{DS})', VGSopt, VDSopt);
plotHeat(Vgs_vec, Vds_vec, GMIDg, 'g_m/I_D [1/V]',   'Heatmap: g_m/I_D(V_{GS},V_{DS})', VGSopt, VDSopt);
plotHeat(Vgs_vec, Vds_vec, AVg,   'A_v=g_m/g_{ds}',  'Heatmap: A_v(V_{GS},V_{DS})', VGSopt, VDSopt);
FEASg = (GMg >= specs.gm_min) & ...
        (GMIDg >= specs.eta_gmId) & ...
        (AVg >= specs.Av_min) & ...
        (VDSg >= specs.Vds_min);

figure;
imagesc(Vgs_vec, Vds_vec, FEASg); 
axis xy; 
grid on;
colormap(parula); 
colorbar;
xlabel('V_{GS} [V]'); 
ylabel('V_{DS} [V]');
ttl = sprintf(['Feasibility Region (1=feasible,0=infeasible)\n' ...
               'g_m >= %.1f mS,  g_m/I_D >= %.1f 1/V,  A_v >= %.1f,  V_DS >= %.1f V'], ...
               1e3*specs.gm_min, specs.eta_gmId, specs.Av_min, specs.Vds_min);
t = title(ttl);
set(t,'Interpreter','none'); 

hold on;
plot(VGSopt, VDSopt, 'rx', 'MarkerSize', 12, 'LineWidth', 2);


VGS_h = 0.70;
VDS_h = 0.60;
ver_opt = verify_bias_point(VGSopt, VDSopt, p, specs);
ver_h   = verify_bias_point(VGS_h,  VDS_h,  p, specs);

ReportTbl = table( ...
    [ver_opt.Vgs; ver_h.Vgs], ...
    [ver_opt.Vds; ver_h.Vds], ...
    [ver_opt.Id;  ver_h.Id], ...
    [ver_opt.gm;  ver_h.gm], ...
    [ver_opt.gds; ver_h.gds], ...
    [ver_opt.gmId;ver_h.gmId], ...
    [ver_opt.Av;  ver_h.Av], ...
    [ver_opt.P;   ver_h.P], ...
    'VariableNames', {'VGS_V','VDS_V','ID_A','gm_S','gds_S','gmId_1perV','Av','P_W'}, ...
    'RowNames', {'Optimized','Heuristic'} );

% sensitivity
gm_sweep = linspace(0.5e-3, 2.0e-3, 10);  % 0.5mS to 2mS
P_sweep  = NaN(size(gm_sweep));
Id_sweep = NaN(size(gm_sweep));
ok       = false(size(gm_sweep));

spec0 = specs;

for k = 1:numel(gm_sweep)
    specs_k = spec0;
    specs_k.gm_min = gm_sweep(k);

    try
        sol_k = solve_bias_lp_cvx(T, specs_k);
        if strcmp(sol_k.cvx_status, "Solved")
            ver_k = verify_bias_point(sol_k.VGS_pick, sol_k.VDS_pick, p, specs_k);
            ok(k) = true;
            P_sweep(k)  = ver_k.P;
            Id_sweep(k) = ver_k.Id;
         end
    catch
      ok(k) = false;
    end
end

figure; grid on; hold on;
plot(1e3*gm_sweep, 1e3*P_sweep, '-o', 'LineWidth', 1.5);
xlabel('g_{m,min} [mS]');
ylabel('Min Power [mW]');
title('Sensitivity: Min Power vs g_{m,min} (eta fixed)');
plot(1e3*gm_sweep(~ok), zeros(sum(~ok),1), 'rx', 'MarkerSize', 10, 'LineWidth', 2);

figure; grid on; hold on;
plot(1e3*gm_sweep, 1e3*Id_sweep, '-o', 'LineWidth', 1.5);
xlabel('g_{m,min} [mS]');
ylabel('Min I_D [mA]');
title('Sensitivity: Min I_D vs g_{m,min} (eta fixed)');
plot(1e3*gm_sweep(~ok), zeros(sum(~ok),1), 'rx', 'MarkerSize', 10, 'LineWidth', 2);

function plotHeat(Vgs_vec, Vds_vec, Z, cbarLabel, ttl, VGSopt, VDSopt)
    figure;
    imagesc(Vgs_vec, Vds_vec, Z); axis xy; grid on;
    cb = colorbar; ylabel(cb, cbarLabel);
    xlabel('V_{GS} [V]'); ylabel('V_{DS} [V]');
    title(ttl);
    hold on;
    plot(VGSopt, VDSopt, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
end

outDir = "report_figures";
if ~exist(outDir, "dir")
     mkdir(outDir);
end

figs = findall(0, 'Type', 'figure');

for k = 1:numel(figs)
    f = figs(k);
    nm = string(get(f,'Name'));
    if strlength(nm) == 0
       nm = "Fig_" + string(f.Number);
    end
    nm = regexprep(nm, '[^\w\- ]', '');
    nm = strrep(nm, ' ', '_');
    savefig(f, fullfile(outDir, nm + ".fig"));
    exportgraphics(f, fullfile(outDir, nm + ".png"), 'Resolution', 300);
end

fprintf("Saved %d figures to folder: %s\n", numel(figs), outDir);
