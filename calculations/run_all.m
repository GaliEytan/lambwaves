function run_all
%RUN_ALL Reproduce the first numerical study for the manuscript.
% Uses only core MATLAB functions and writes all results under results/.

root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root, 'results');
if ~exist(out, 'dir'), mkdir(out); end

mat.name = 'Aluminum (representative isotropic values)';
mat.rho = 2700;            % kg/m^3
mat.cL = 6320;             % m/s
mat.cT = 3130;             % m/s
mat.d = 1.0e-3;            % total plate thickness, m
mat.h = mat.d/2;

opts.K = linspace(0.04, 8.0, 260).';       % dimensionless kh
opts.Wmax = 12.0;                          % omega*h/cT
opts.scanPoints = 7000;
opts.nBranches = 3;

fprintf('Solving Rayleigh-Lamb equations...\n');
[A, Ares] = solve_family('A', mat, opts);
[S, Sres] = solve_family('S', mat, opts);

% Dimensional quantities and derivatives.
A = add_kinematics(A, opts.K, mat);
S = add_kinematics(S, opts.K, mat);

% Low-frequency validation.
E = mat.rho*mat.cT^2*(3*mat.cL^2-4*mat.cT^2)/(mat.cL^2-mat.cT^2);
nu = (mat.cL^2-2*mat.cT^2)/(2*(mat.cL^2-mat.cT^2));
cPlate = sqrt(E/(mat.rho*(1-nu^2)));
Dplate = E*mat.d^3/(12*(1-nu^2));
kSmall = opts.K(1:15)/mat.h;
S0expected = cPlate*kSmall;
A0expected = sqrt(Dplate/(mat.rho*mat.d))*kSmall.^2;
S0err = median(abs(S.omega(1:15,1)-S0expected)./S0expected);
A0err = median(abs(A.omega(1:15,1)-A0expected)./A0expected);

% ZGV candidates are local frequency minima away from endpoints.
zgvA = find_zgv(A, opts.K, mat, 'A');
zgvS = find_zgv(S, opts.K, mat, 'S');

% Causal reduced-order calculation on the fundamental branches.
fprintf('Evaluating causal modal integrals...\n');
transient = causal_study(A, S, opts.K, mat);

% Tables and plots.
write_dispersion_csv(fullfile(out, 'dispersion.csv'), opts.K, A, S, mat);
plot_dispersion(fullfile(out, 'dispersion.png'), A, S, mat);
plot_group_velocity(fullfile(out, 'group_velocity.png'), A, S, mat);
plot_transient(fullfile(out, 'causal_buildup.png'), transient);
write_summary(fullfile(out, 'validation_summary.txt'), mat, opts, Ares, Sres, ...
    S0err, A0err, cPlate, zgvA, zgvS, transient);

save(fullfile(out, 'calculation_data.mat'), 'mat', 'opts', 'A', 'S', ...
    'Ares', 'Sres', 'zgvA', 'zgvS', 'transient');
fprintf('Finished. Results written to %s\n', out);
end

function [branches, maxResidual] = solve_family(family, mat, opts)
K = opts.K;
Wgrid = linspace(1e-4, opts.Wmax, opts.scanPoints);
W = nan(numel(K), opts.nBranches);
residuals = nan(size(W));
for i = 1:numel(K)
    vals = arrayfun(@(w) characteristic(w, K(i), mat.cT/mat.cL, family), Wgrid);
    roots = [];
    for j = 1:numel(Wgrid)-1
        if near_bulk_threshold(Wgrid(j), K(i), mat.cT/mat.cL), continue; end
        if isfinite(vals(j)) && isfinite(vals(j+1)) && vals(j)*vals(j+1) < 0
            try
                r = fzero(@(w) characteristic(w, K(i), mat.cT/mat.cL, family), ...
                    [Wgrid(j), Wgrid(j+1)]);
                if ~near_bulk_threshold(r, K(i), mat.cT/mat.cL) && ...
                        (isempty(roots) || min(abs(roots-r)) > 2e-4)
                    roots(end+1) = r; %#ok<AGROW>
                end
            catch
            end
        end
    end
    roots = sort(roots);
    take = min(opts.nBranches, numel(roots));
    if take > 0
        W(i,1:take) = roots(1:take);
        for b = 1:take
            residuals(i,b) = normalized_residual(roots(b), K(i), mat.cT/mat.cL, family);
        end
    end
end
branches.W = W;
branches.omega = W*mat.cT/mat.h;
maxResidual = max(residuals(:), [], 'omitnan');
end

function tf = near_bulk_threshold(W, K, r)
tf = abs(W-K) < 5e-3 || abs(W-K/r) < 5e-3;
end

function y = characteristic(W, K, r, family)
p = sqrt(complex((r*W)^2-K^2));
q = sqrt(complex(W^2-K^2));
if family == 'S'
    detv = (q^2-K^2)^2*cos(p)*sin(q) + 4*K^2*p*q*sin(p)*cos(q);
else
    detv = 4*K^2*p*q*cos(p)*sin(q) + (q^2-K^2)^2*sin(p)*cos(q);
end
if W < K
    y = imag(detv);
else
    y = real(detv);
end
scale = 1 + abs((q^2-K^2)^2) + abs(4*K^2*p*q);
y = y/scale;
end

function r = normalized_residual(W, K, ratio, family)
r = abs(characteristic(W, K, ratio, family));
end

function B = add_kinematics(B, K, mat)
k = K/mat.h;
B.k = k;
for b = 1:size(B.omega,2)
    good = isfinite(B.omega(:,b));
    B.vg(:,b) = nan(size(k));
    B.curvature(:,b) = nan(size(k));
    if nnz(good) > 4
        B.vg(good,b) = gradient(B.omega(good,b), k(good));
        B.curvature(good,b) = gradient(B.vg(good,b), k(good));
    end
end
end

function z = find_zgv(B, K, mat, family)
z = struct('mode', {}, 'K', {}, 'fd_MHz_mm', {}, 'frequency_MHz', {});
for b = 2:size(B.W,2) % exclude fundamental acoustic branch
    w = B.W(:,b);
    for i = 3:numel(K)-2
        if all(isfinite(w(i-1:i+1))) && w(i) < w(i-1) && w(i) < w(i+1)
            n = numel(z)+1;
            z(n).mode = sprintf('%s%d', family, b-1);
            z(n).K = K(i);
            z(n).fd_MHz_mm = B.omega(i,b)/(2*pi)*mat.d/1000;
            z(n).frequency_MHz = B.omega(i,b)/(2*pi)/1e6;
        end
    end
end
end

function R = causal_study(A, S, K, mat)
% Normalized modal coupling: Gaussian line source and scalar modal response.
targetK = 1.0;
[~, idx] = min(abs(K-targetK));
Omega = A.omega(idx,1);
w0 = 1.0e-3;
x = 20e-3;
gamma = 0.01*Omega;
k = A.k;
mask = isfinite(A.omega(:,1)) & isfinite(S.omega(:,1));
k = k(mask); Kuse = K(mask);
wa = A.omega(mask,1); ws = S.omega(mask,1);
source = exp(-(k*w0).^2/4);
dk = mean(diff(k));
weights = source.*cos(k*x)*dk/pi;

T = 2*pi/Omega;
t = linspace(0, 250*T, 2501).';
[ua, usa, utra] = modal_integral(wa, weights, Omega, gamma, t);
[us, uss, utrs] = modal_integral(ws, weights, Omega, gamma, t);
ratioA = abs(utra)/max(abs(usa), eps);
futureMax = flipud(cummax(flipud(ratioA)));
j = find(futureMax < 0.05, 1, 'first');
if isempty(j), tss = NaN; else, tss = t(j); end

vgTarget = A.vg(idx,1);
tArrival = x/abs(vgTarget);
N95 = ceil(tss/T);
purityA = abs(ua).^2./(abs(ua).^2+abs(us).^2+eps);

R.t = t; R.T = T; R.frequency = Omega/(2*pi); R.Omega = Omega;
R.uA = ua; R.uS = us; R.transientA = utra; R.ssA = usa;
R.ratioA = ratioA; R.purityA = purityA;
R.tss = tss; R.tArrival = tArrival; R.N95 = N95;
R.x = x; R.w0 = w0; R.gamma = gamma;
R.targetK = Kuse(find(mask(idx:end),1,'first')); %#ok<FNDSB>
end

function [u, ussTime, utr] = modal_integral(wn, weights, Omega, gamma, t)
H = weights./(wn.^2-Omega^2-2i*gamma*Omega);
wd = sqrt(max(wn.^2-gamma^2, 0));
sp = -gamma+1i*wd;
sm = -gamma-1i*wd;
Cp = H.*(1i*Omega+sm)./(sp-sm);
Cm = -H-Cp;
ssPhasor = sum(H);
ussTime = ssPhasor*exp(-1i*Omega*t);
utr = zeros(size(t));
for j = 1:numel(wn)
    utr = utr + Cp(j)*exp(sp(j)*t) + Cm(j)*exp(sm(j)*t);
end
u = ussTime+utr;
end

function write_dispersion_csv(path, K, A, S, mat)
fid = fopen(path, 'w');
fprintf(fid, 'K');
for b=1:3, fprintf(fid, ',A%d_fd_MHz_mm,A%d_vg_m_s',b-1,b-1); end
for b=1:3, fprintf(fid, ',S%d_fd_MHz_mm,S%d_vg_m_s',b-1,b-1); end
fprintf(fid, '\n');
for i=1:numel(K)
    fprintf(fid, '%.9g', K(i));
    for b=1:3, fprintf(fid, ',%.9g,%.9g', A.omega(i,b)/(2*pi)*mat.d/1000, A.vg(i,b)); end
    for b=1:3, fprintf(fid, ',%.9g,%.9g', S.omega(i,b)/(2*pi)*mat.d/1000, S.vg(i,b)); end
    fprintf(fid, '\n');
end
fclose(fid);
end

function plot_dispersion(path, A, S, mat)
f = figure('Visible','off','Color','w','Position',[100 100 900 560]); hold on;
colors = lines(3);
for b=1:3
    plot(A.k*mat.d, A.omega(:,b)/(2*pi)*mat.d/1000, '--', 'Color',colors(b,:), 'LineWidth',1.5);
    plot(S.k*mat.d, S.omega(:,b)/(2*pi)*mat.d/1000, '-', 'Color',colors(b,:), 'LineWidth',1.5);
end
xlabel('k d'); ylabel('f d (MHz mm)'); grid on; box on;
legend('A_0','S_0','A_1','S_1','A_2','S_2','Location','northwest');
title('Rayleigh-Lamb dispersion: 1 mm isotropic aluminum plate');
exportgraphics(f,path,'Resolution',180); close(f);
end

function plot_group_velocity(path, A, S, mat)
f = figure('Visible','off','Color','w','Position',[100 100 900 560]); hold on;
plot(A.omega(:,1)/(2*pi)*mat.d/1000,A.vg(:,1),'-','LineWidth',1.6);
plot(S.omega(:,1)/(2*pi)*mat.d/1000,S.vg(:,1),'-','LineWidth',1.6);
yline(0,'k:'); xlabel('f d (MHz mm)'); ylabel('group velocity (m/s)');
legend('A_0','S_0','Location','best'); grid on; box on;
title('Fundamental-mode group velocity'); exportgraphics(f,path,'Resolution',180); close(f);
end

function plot_transient(path, R)
f = figure('Visible','off','Color','w','Position',[100 100 900 700]);
tcycles=R.t/R.T;
subplot(2,1,1); plot(tcycles,abs(R.uA)/max(abs(R.ssA)),'LineWidth',1.3); hold on;
yline(1,'k:'); xlabel('time (drive cycles)'); ylabel('|u_{A0}|/|u_{ss}|'); grid on;
title('Normalized causal buildup (scalar modal coupling)');
subplot(2,1,2); plot(tcycles,R.purityA,'LineWidth',1.3);
xlabel('time (drive cycles)'); ylabel('A_0 modal purity'); ylim([0 1]); grid on;
exportgraphics(f,path,'Resolution',180); close(f);
end

function write_summary(path, mat, opts, Ares, Sres, S0err, A0err, cPlate, zA, zS, R)
fid=fopen(path,'w');
fprintf(fid,'CALCULATION AND VALIDATION SUMMARY\n');
fprintf(fid,'Material: %s\n',mat.name);
fprintf(fid,'rho=%.1f kg/m^3, cL=%.1f m/s, cT=%.1f m/s, d=%.6g m\n',mat.rho,mat.cL,mat.cT,mat.d);
fprintf(fid,'Grid: %d K values, %d frequency scan points, Wmax=%.2f\n\n',numel(opts.K),opts.scanPoints,opts.Wmax);
fprintf(fid,'Maximum normalized characteristic residual: A %.3e, S %.3e\n',Ares,Sres);
fprintf(fid,'Low-frequency median relative error: A0 thin-plate %.3f%%, S0 plate-speed %.3f%%\n',100*A0err,100*S0err);
fprintf(fid,'Expected low-frequency S0 plate speed: %.3f m/s\n\n',cPlate);
fprintf(fid,'Detected higher-mode local minima (grid estimates):\n');
allz=[zA,zS];
if isempty(allz), fprintf(fid,'  none on sampled branches/range\n'); end
for i=1:numel(allz)
    fprintf(fid,'  %s: K=%.4f, f*d=%.5f MHz mm, f=%.5f MHz\n',allz(i).mode,allz(i).K,allz(i).fd_MHz_mm,allz(i).frequency_MHz);
end
fprintf(fid,'\nNormalized causal study (not an absolute thermoelastic amplitude):\n');
fprintf(fid,'  drive f=%.6f MHz, x=%.3f mm, Gaussian width=%.3f mm, gamma/Omega=%.4f\n',R.frequency/1e6,R.x*1e3,R.w0*1e3,R.gamma/R.Omega);
fprintf(fid,'  group-delay estimate x/|vg|=%.6g us\n',R.tArrival*1e6);
fprintf(fid,'  5%% envelope settling time=%.6g us (%.1f cycles)\n',R.tss*1e6,R.tss/R.T);
fprintf(fid,'  N95 proxy (cycles to remain within 5%%)=%.0f\n',R.N95);
fprintf(fid,'\nValidation scope:\n');
fprintf(fid,'  Dispersion results solve the full isotropic Rayleigh-Lamb equations.\n');
fprintf(fid,'  The causal calculation uses unit scalar modal overlap and constant damping; it validates the computational framework, not an experimental absolute amplitude.\n');
fclose(fid);
end
