%% RUN_PAPER_EXPERIMENTS  Reproduce the experiments and figures of the paper.
% Experiments of Section "Computational Cost and Accuracy of RSH Rotation".
%
% Usage (from MATLAB, with the repository root as current folder):
%   run('experiments/run_paper_experiments.m')
% The script adds src/ to the path, writes the five figures of the paper to
% figures/*.pdf and the numerical results to
% results/rsh_rot_experiments_results.mat (both inside the repository
% folder, created if missing). Set cfg.saveFigs = false to
% only display the figures. A full run takes from several minutes to about an
% hour depending on the machine (l up to 64); reduce cfg.lList for a quick run.
% Random draws are reproducible (rng(3)); timings depend on the hardware.
%
% Methods
%   'Ivanic'  Ivanic-Ruedenberg recursion (reference). Builds all bands 1..l.
%   'Lie-exp' M_l = expm(theta * n.L^l)                    (generic matrix exponential)
%   'Lie-eig' M_l = Re(V diag(exp(lambda)) V^-1), [V,lambda] = eig(theta * n.L^l)
%   'Lie-AD'  anti-diagonal factorization M_l = F_l * M_lz(theta) * F_l'
%             (proposed method: rsh_rot_lie_eig / rsh_rot_to_z). Builds all bands 1..l.
%
% Requires on the path
%   rsh_lie_generator(l)   -> [Lx, Ly, Lz], each (2l+1)x(2l+1), band l only
%   rsh_rot_lie_eig(lmax, x, gens) -> [Ml, M] (proposed method, x = Rodrigues vector)
%   rsh_rot_to_z(l, xhat, G) -> F_l, (2l+1)x(2l+1), used by rsh_rot_lie_eig
%   (versions that accept precomputed generators, so that Lie-AD, like Lie-exp and
%   Lie-eig, does not rebuild them inside the timed loops)
% No toolbox is needed (axis-angle conversions are implemented below).
%
% Experiments
%   0. Convention checks: equivariance y_l(R u) = M_l y_l(u) for every method, which
%      fixes the exponent sign, the real-SH phase convention and the Ivanic variant
%      independently of any reference; self-check of F_l.
%   1. Time per rotation vs l, (a) for ALL bands 1..l, (b) for band l only.
%   2. Discrepancy from Ivanic, Lie-exp vs Lie-AD, and orthogonality residual.
%   3. Rotation of von Mises-Fisher distributions with closed-form coefficients.
%   4. Many angles about a fixed axis: reuse of F_l (coefficients and correlation C(theta)).
%   5. Closed-form gradient/Hessian vs finite differences; Newton alignment on SO(3).

%% 0. Configuration
clear; clc; rng(3);
% Repository root. mfilename is empty when the file is evaluated section by
% section (Run Section / Ctrl+Enter), so fall back to WHICH and to the
% current folder; output paths below are then always absolute.
repoCand = {fileparts(fileparts(mfilename('fullpath'))), ...
            fileparts(fileparts(which('run_paper_experiments'))), pwd, fileparts(pwd)};
repoDir = '';
for iCand = 1:numel(repoCand)
    if ~isempty(repoCand{iCand}) && exist(fullfile(repoCand{iCand}, 'src', 'rsh_lie_generator.m'), 'file') == 2
        repoDir = repoCand{iCand}; break;
    end
end
if isempty(repoDir)
    error('Cannot locate the repository root: cd to the repository (or to experiments/) and run again.');
end
clear repoCand iCand
addpath(genpath(fullfile(repoDir, 'src')));
cfg.lList       = [1 2 3 4 6 8 12 16 24 32 48 64];  % bands to test
cfg.nRotTime    = 100;   % rotations per timing loop for small l (time reported per rotation)
cfg.nRotTimeMin = 10;    % fewer rotations for large l, to keep the run time reasonable
cfg.nRotAcc     = 50;    % random rotations per l in Experiment 2
cfg.nTrialVMF   = 20;    % trials per configuration in Experiment 3
cfg.includeConversion = true;  % true: every method receives the 3x3 matrix R, so the Lie
                               % methods also pay for the axis-angle extraction
cfg.tolConv     = 1e-8;         % tolerance of the convention checks
cfg.saveFigs    = true;        % save every figure as a vector PDF
cfg.figDir      = fullfile(repoDir, 'figures');   % output folder of the PDFs
cfg.resDir      = fullfile(repoDir, 'results');   % output folder of the .mat results
cfg.vmf = struct('name',  {'\kappa=5', '\kappa=50', '\kappa=500', 'mixture'}, ...
                 'w',     {1, 1, 1, [0.5 0.3 0.2]}, ...
                 'kappa', {5, 50, 500, [20 80 300]});

% --- Additional experiments (4-5): properties specific to the Lie-algebraic formulation
cfg.mix        = struct('w', [0.5 0.3 0.2], 'kappa', [20 80 300]);  % test distribution
cfg.fixAxisL   = [8 16 32 64];   % Test 4: bands 0..L
cfg.fixAxisK   = 360;            % Test 4: angles per sweep (1 degree step)
cfg.fixAxisTrials = 5;           % Test 4: random axes / distributions per L (accuracy)
cfg.gradL      = 16;             % Test 5: bands 0..L of the correlation function
cfg.fdSteps    = 10.^(-(1:0.5:8));  % Test 5: finite-difference steps
cfg.newtonTrials = 20;           % Test 5: Newton alignments
cfg.newtonInit = 30*pi/180;      % Test 5: initial rotation error [rad]
cfg.newtonIters = 20;            % Test 5: iterations

nL   = numel(cfg.lList);
lMax = max(cfg.lList);
methods = {'Ivanic', 'Lie-exp', 'Lie-eig', 'Lie-AD'};

% Generators of all bands 1..lMax, computed once and NOT included in the timings
% (they depend only on l, not on the rotation).
gens = cell(1, lMax);
for l = 1:lMax
    [gx, gy, gz] = rsh_lie_generator(l);
    assert(isequal(size(gx), [2*l+1, 2*l+1]), 'rsh_lie_generator(%d) has wrong size', l);
    gens{l} = struct('Lx', gx, 'Ly', gy, 'Lz', gz);
end

% ======================================================================
%% 0. CONVENTION CHECKS (equivariance, no reference needed)
% ======================================================================
% For each method the residual max ||y_l(R u) - M_l y_l(u)|| over random points u and
% rotations R is computed for l = 2, 3 and for each possible convention:
%   cs    : real SH evaluated without (0) or with (1) the (-1)^m phase
%   sgn   : exponent sign of the Lie methods, M = expm(sgn*theta*n.L)
%   sgnAD : sign of the Rodrigues vector passed to rsh_rot_lie_eig
%   ivMode: Ivanic variant, 1: M(R), 2: M(R'), 3: D M(R) D, 4: D M(R') D, D = diag((-1)^m)
% The generators define the convention of the paper, so cs and sgn are chosen from
% Lie-exp; the other methods are then checked against the same evaluator.
nU = 30; nR = 10; lc = [2 3];
Uc = rand_unit(nU); Rc = cell(1,nR); Xc = zeros(3,nR);
for k = 1:nR
    [ax, ang] = rand_axang(); Rc{k} = axang_to_rotm(ax, ang); Xc(:,k) = ax*ang;
end
resExp = zeros(2,2);                     % (cs+1, sign index) with sign index 1:+, 2:-
for cs = 0:1
    for is = 1:2
        s = 3 - 2*is;
        r = 0;
        for l = lc
            for k = 1:nR
                A = lie_alg(gens{l}, Xc(:,k));
                r = max(r, equiv_res(expm(s*A), l, Rc{k}, Uc, cs));
            end
        end
        resExp(cs+1, is) = r;
    end
end
[resLie, idx] = min(resExp(:));
[ics, isg] = ind2sub([2 2], idx);
csConv = ics - 1; sgn = 3 - 2*isg;
fprintf('Lie-exp equivariance: cs=0 (+/-): %.1e %.1e | cs=1 (+/-): %.1e %.1e -> cs=%d, sgn=%+d\n', ...
    resExp(1,1), resExp(1,2), resExp(2,1), resExp(2,2), csConv, sgn);
if resLie > cfg.tolConv
    warning(['Lie-exp is not equivariant for any convention: the generators and the ' ...
             'real-SH evaluator rsh_eval_ref do not use the same basis.']);
end

resAD = zeros(1,2);
for is = 1:2
    s = 3 - 2*is; r = 0;
    for l = lc
        for k = 1:nR
            Mad = rsh_rot_lie_eig(l, s*Xc(:,k), gens);
            r = max(r, equiv_res(Mad, l, Rc{k}, Uc, csConv));
        end
    end
    resAD(is) = r;
end
[resADmin, isg] = min(resAD); sgnAD = 3 - 2*isg;
fprintf('Lie-AD  equivariance: +x %.1e | -x %.1e -> sgnAD=%+d\n', resAD(1), resAD(2), sgnAD);
if resADmin > cfg.tolConv
    warning('Lie-AD (rsh_rot_lie_eig) is not equivariant in the convention of the generators.');
end

resIv = zeros(1,4);
for mode = 1:4
    r = 0;
    for l = lc
        for k = 1:nR
            r = max(r, equiv_res(ivanic_band(l, Rc{k}, mode), l, Rc{k}, Uc, csConv));
        end
    end
    resIv(mode) = r;
end
[resIvMin, ivMode] = min(resIv);
fprintf('Ivanic  equivariance: M(R) %.1e | M(R'') %.1e | DM(R)D %.1e | DM(R'')D %.1e -> mode %d\n', ...
    resIv, ivMode);
if resIvMin > cfg.tolConv
    warning('Ivanic is not equivariant in any variant: check the recursion formulas.');
end

% Self-check of the similarity matrix on all tested bands:
% F_l orthogonal and F_l' L_{xhat} F_l = L_z.
errF = zeros(1, nL); errFo = errF;
for il = 1:nL
    l = cfg.lList(il);
    for k = 1:nR
        xh = Xc(:,k)/norm(Xc(:,k));
        F  = rsh_rot_to_z(l, sgnAD*xh, gens{l});
        errFo(il) = max(errFo(il), norm(F'*F - eye(2*l+1), 'fro'));
        errF(il)  = max(errF(il),  norm(F'*lie_alg(gens{l}, sgnAD*xh)*F - gens{l}.Lz, 'fro'));
    end
    fprintf('F_l self-check l=%2d: max ||F''F-I|| = %.1e, max ||F''L F - Lz|| = %.1e\n', l, errFo(il), errF(il));
end

% ======================================================================
%% TEST 1: computation time vs l
% ======================================================================
% (a) T_all : time to obtain ALL bands 1..l (rotating a full expansion). Ivanic and
%             Lie-AD produce all bands by construction; Lie-exp and Lie-eig loop over bands.
% (b) T_band: time to obtain band l only. Not available for Ivanic, which needs the
%             lower bands; Lie-AD uses ad_band (same algorithm, band l only).
% The same rotations are given to all methods. Times are per rotation.
T_all  = nan(nL, 4);
T_band = nan(nL, 4);
for il = 1:nL
    l = cfg.lList(il);
    nRot = max(cfg.nRotTimeMin, round(cfg.nRotTime * min(1, (16/l)^2)));
    Rl = cell(nRot,1); Xl = zeros(3, nRot);
    for k = 1:nRot
        [ax, ang] = rand_axang(); Rl{k} = axang_to_rotm(ax, ang); Xl(:,k) = ax*ang;
    end
    conv = cfg.includeConversion;
    T_all(il,1)  = timeit(@() run_ivanic(l, Rl, ivMode))                     / nRot;
    T_all(il,2)  = timeit(@() run_lie(gens(1:l), Rl, Xl, 'exp', sgn, conv))   / nRot;
    T_all(il,3)  = timeit(@() run_lie(gens(1:l), Rl, Xl, 'eig', sgn, conv))   / nRot;
    T_all(il,4)  = timeit(@() run_ad(l, Rl, Xl, sgnAD, conv, true, gens))     / nRot;
    T_band(il,2) = timeit(@() run_lie(gens(l), Rl, Xl, 'exp', sgn, conv))     / nRot;
    T_band(il,3) = timeit(@() run_lie(gens(l), Rl, Xl, 'eig', sgn, conv))     / nRot;
    T_band(il,4) = timeit(@() run_ad(l, Rl, Xl, sgnAD, conv, false, gens))    / nRot;
    fprintf(['l=%2d | all bands [s]: Ivanic %.2e  Lie-exp %.2e  Lie-eig %.2e  Lie-AD %.2e' ...
             ' | band l: Lie-exp %.2e  Lie-eig %.2e  Lie-AD %.2e\n'], l, T_all(il,:), T_band(il,2:4));
end

c = lines(4); mk = {'o','s','^','d'};
fig = new_fig(12, 9);
for j = 1:4
    loglog(cfg.lList, T_all(:,j), ['-' mk{j}], 'Color', c(j,:), 'LineWidth', 1.2, ...
        'DisplayName', [methods{j} ' (bands 1..l)']); hold on
end
for j = 2:4
    loglog(cfg.lList, T_band(:,j), ['--' mk{j}], 'Color', c(j,:), 'LineWidth', 1.0, ...
        'DisplayName', [methods{j} ' (band l)']);
end
set(gca, 'XScale', 'log', 'YScale', 'log'); grid on
xlabel('band l'); ylabel('time per rotation [s]'); legend('show', 'Location', 'northwest');
title('Computation time vs l');
save_fig(fig, cfg, 'rsh_time');


% ======================================================================
%% TEST 2: discrepancy from Ivanic and orthogonality vs l
% ======================================================================
% E(:,:,1..3): ||M_X - M_Ivanic||_F for X = Lie-exp, Lie-eig, Lie-AD
% E_expAD    : ||M_Lie-exp - M_Lie-AD||_F (two Lie methods, independent of the recursion)
% O(:,:,1..4): ||M_X' M_X - I||_F for all four methods (reference-free accuracy check)
E = nan(nL, cfg.nRotAcc, 3); E_expAD = nan(nL, cfg.nRotAcc); O = nan(nL, cfg.nRotAcc, 4);
for il = 1:nL
    l = cfg.lList(il); I = eye(2*l+1);
    for k = 1:cfg.nRotAcc
        [ax, ang] = rand_axang(); R = axang_to_rotm(ax, ang); x = ax*ang;
        Ms = all_methods(l, R, x, gens, sgn, sgnAD, ivMode);
        for j = 1:3, E(il,k,j) = norm(Ms{j+1} - Ms{1}, 'fro'); end
        E_expAD(il,k) = norm(Ms{2} - Ms{4}, 'fro');
        for j = 1:4, O(il,k,j) = norm(Ms{j}'*Ms{j} - I, 'fro'); end
    end
    fprintf(['l=%2d | median ||M-Ivanic||_F: Lie-exp %.1e  Lie-eig %.1e  Lie-AD %.1e' ...
             ' | Lie-exp vs Lie-AD %.1e | median ||M''M-I||_F: %.1e %.1e %.1e %.1e\n'], ...
        l, median(E(il,:,1)), median(E(il,:,2)), median(E(il,:,3)), median(E_expAD(il,:)), ...
        median(O(il,:,1)), median(O(il,:,2)), median(O(il,:,3)), median(O(il,:,4)));
end

fig = new_fig(18, 8);
subplot(1,2,1);
for j = 1:3
    loglog(cfg.lList, median(E(:,:,j),2) + eps, ['-' mk{j+1}], 'Color', c(j+1,:), ...
        'DisplayName', [methods{j+1} ' vs Ivanic (median)']); hold on
    loglog(cfg.lList, max(E(:,:,j),[],2) + eps, ['--' mk{j+1}], 'Color', c(j+1,:), ...
        'DisplayName', [methods{j+1} ' vs Ivanic (max)']);
end
loglog(cfg.lList, median(E_expAD,2) + eps, '-x', 'Color', [0 0 0], 'DisplayName', 'Lie-exp vs Lie-AD (median)');
set(gca, 'XScale', 'log', 'YScale', 'log'); grid on
xlabel('band l'); ylabel('||M - M_{ref}||_F'); legend('show', 'Location', 'northwest');
title('Discrepancy vs l');
subplot(1,2,2);
for j = 1:4
    loglog(cfg.lList, median(O(:,:,j),2) + eps, ['-' mk{j}], 'Color', c(j,:), ...
        'DisplayName', methods{j}); hold on
end
set(gca, 'XScale', 'log', 'YScale', 'log'); grid on
xlabel('band l'); ylabel('||M^T M - I||_F (median)'); legend('show', 'Location', 'northwest');
title('Orthogonality residual vs l');
save_fig(fig, cfg, 'rsh_accuracy');


% ======================================================================
%% TEST 3: rotation of von Mises-Fisher distributions
% ======================================================================
% s_l = sum_k w_k g_l(kappa_k) y_l(mu_k),  d_l = sum_k w_k g_l(kappa_k) y_l(R mu_k),
% err_X(l) = ||d_l - M_X s_l||_2. The same rotations and matrices are used for all
% configurations. Absolute and relative errors are stored.
nC = numel(cfg.vmf);
errA = nan(nC, nL, cfg.nTrialVMF, 4);   % absolute error
dnrm = nan(nC, nL, cfg.nTrialVMF);      % ||d_l||_2
for t = 1:cfg.nTrialVMF
    [ax, ang] = rand_axang(); R = axang_to_rotm(ax, ang); x = ax*ang;
    Mi_all = ivanic_all(lMax, R, ivMode);
    [~, Mad_all] = rsh_rot_lie_eig(lMax, sgnAD*x, gens);
    mus = cell(1,nC);
    for ic = 1:nC, mus{ic} = rand_unit(numel(cfg.vmf(ic).kappa)); end
    for il = 1:nL
        l = cfg.lList(il);
        A  = lie_alg(gens{l}, sgn*x);
        Ms = {Mi_all{l}, expm(A), lie_eig(A), Mad_all{l+1}};
        for ic = 1:nC
            v  = cfg.vmf(ic);
            gl = v.w(:)' .* vmf_g(l, v.kappa(:)');               % 1 x K
            s  = rsh_eval_ref(l, mus{ic},   csConv) * gl';
            d  = rsh_eval_ref(l, R*mus{ic}, csConv) * gl';
            dnrm(ic, il, t) = norm(d);
            for j = 1:4, errA(ic, il, t, j) = norm(d - Ms{j}*s); end
        end
    end
end
errR = errA ./ repmat(dnrm, [1 1 1 4]);
errR(repmat(dnrm, [1 1 1 4]) < realmin) = NaN;

% Check of the coefficients: sum_l (2l+1) g_l(kappa)/(4 pi) = f(mu; mu, kappa)
for kap = [5 50 500]
    lS = 0:max(200, ceil(10*sqrt(kap)));
    fS = sum((2*lS+1) .* vmf_g(lS, kap)) / (4*pi);
    f0 = kap / (2*pi*(1 - exp(-2*kap)));                      % kappa e^k / (4 pi sinh k)
    fprintf('vMF check kappa=%3d: sum over %d bands %.10e, f(mu) %.10e, rel. diff %.1e\n', ...
        kap, numel(lS), fS, f0, abs(fS-f0)/f0);
end
for ic = 1:nC
    for il = 1:nL
        fprintf('%-10s l=%2d | ||d||=%.1e | median rel. err: Ivanic %.1e  Lie-exp %.1e  Lie-eig %.1e  Lie-AD %.1e\n', ...
            cfg.vmf(ic).name, cfg.lList(il), median(dnrm(ic,il,:)), ...
            median(errR(ic,il,:,1)), median(errR(ic,il,:,2)), median(errR(ic,il,:,3)), median(errR(ic,il,:,4)));
    end
end

fig = new_fig(18, 14);
for ic = 1:nC
    subplot(2,2,ic);
    for j = 1:4
        e = squeeze(errA(ic,:,:,j));
        loglog(cfg.lList, median(e,2) + eps*1e-3, ['-' mk{j}], 'Color', c(j,:), ...
            'DisplayName', [methods{j} ' (median)']); hold on
        loglog(cfg.lList, max(e,[],2) + eps*1e-3, ['--' mk{j}], 'Color', c(j,:), ...
            'HandleVisibility', 'off');
    end
    loglog(cfg.lList, median(squeeze(dnrm(ic,:,:)),2), 'k-', 'LineWidth', 1.5, 'DisplayName', '||d_l||_2');
    set(gca, 'XScale', 'log', 'YScale', 'log'); grid on
    xlabel('band l'); ylabel('error'); title(cfg.vmf(ic).name);
    if ic == 1, legend('show', 'Location', 'southwest'); end
end
save_fig(fig, cfg, 'rsh_vmf');


% ======================================================================
%% TEST 4: many angles about a fixed axis (reuse of F_l)
% ======================================================================
% F_l depends only on the axis. Once it is known, rotating about that axis costs only
% the elementary rotations M_lz(theta), so a sweep over K angles is cheap:
%   (a) rotated coefficients M_l(theta_k) s_l = F_l (M_lz(theta_k) (F_l' s_l)): one dense
%       product per band for all angles together;
%   (b) correlation C(theta) = sum_l d_l' M_l(theta) s_l is a trigonometric polynomial
%       a0 + sum_m (A_m cos m theta + B_m sin m theta), whose coefficients are obtained
%       once; it can then be evaluated at any angle and maximized exactly (Newton).
% Ivanic and Lie-exp must build the matrices for every angle. Lie-eig can also reuse
% its eigendecomposition for a fixed axis (complex arithmetic, non-orthogonal V) and is
% included as the natural competitor. Reference values come from the closed-form vMF
% coefficients (no rotation matrix). Times include the per-axis precomputation.
fa_methods = {'Ivanic', 'Lie-exp', 'Lie-eig (reuse)', 'Lie-AD (reuse)'};
nFA = numel(cfg.fixAxisL); K = cfg.fixAxisK; thK = 2*pi*(0:K-1)/K;
T4coef = nan(nFA, 4); T4corr = nan(nFA, 4); T4pre = nan(nFA, 2);   % pre: Lie-eig, Lie-AD
E4coef = zeros(nFA, 4); E4corr = zeros(nFA, 4); E4ang = nan(nFA, cfg.fixAxisTrials, 5);
for iF = 1:nFA
    L = cfg.fixAxisL(iF);
    for t = 1:cfg.fixAxisTrials
        ax = rand_unit(1); mus = rand_unit(numel(cfg.mix.w)); alpha = 2*pi*rand;
        s = mix_coeffs(L, mus, cfg.mix, eye(3), csConv);
        d = mix_coeffs(L, mus, cfg.mix, axang_to_rotm(ax, alpha), csConv);
        St = mix_coeffs_sweep(L, mus, cfg.mix, ax, thK, csConv);       % truth, cell n_l x K
        Ct = zeros(1, K); nrm = 0;
        for l = 0:L, Ct = Ct + d{l+1}' * St{l+1}; nrm = nrm + norm(d{l+1})*norm(s{l+1}); end
        for j = 1:4
            out = fa_coef(fa_methods{j}, L, s, ax, thK, gens, sgn, sgnAD, ivMode);
            for l = 0:L
                if norm(s{l+1}) > 0
                    E4coef(iF,j) = max(E4coef(iF,j), max(sqrt(sum((out{l+1} - St{l+1}).^2,1))) / norm(s{l+1}));
                end
            end
            C = fa_corr(fa_methods{j}, L, s, d, ax, thK, gens, sgn, sgnAD, ivMode);
            E4corr(iF,j) = max(E4corr(iF,j), max(abs(C - Ct)) / nrm);
            E4ang(iF,t,j) = ang_dist(grid_argmax(C, thK), alpha);      % grid + parabola
        end
        % Lie-AD: exact maximization of the trigonometric polynomial (Newton from grid)
        [a0, A, B] = ad_trig_coeffs(L, s, d, sgnAD*ax, gens);
        E4ang(iF,t,5) = ang_dist(trig_argmax(A, B, grid_argmax(trig_eval(a0,A,B,thK), thK)), alpha);
        if t == 1                                                      % timings on trial 1
            for j = 1:4
                T4coef(iF,j) = timeit(@() fa_coef(fa_methods{j}, L, s, ax, thK, gens, sgn, sgnAD, ivMode));
                T4corr(iF,j) = timeit(@() fa_corr(fa_methods{j}, L, s, d, ax, thK, gens, sgn, sgnAD, ivMode));
            end
            T4pre(iF,1) = timeit(@() eig_all_bands(L, ax, gens, sgn));
            T4pre(iF,2) = timeit(@() F_all_bands(L, sgnAD*ax, gens));
        end
    end
    fprintf(['Test 4, L=%2d, K=%d | coef. time [s]: Ivanic %.2e  Lie-exp %.2e  Lie-eig %.2e  Lie-AD %.2e' ...
             ' | corr. time [s]: %.2e %.2e %.2e %.2e | precomp. F (all bands) %.2e\n'], ...
        L, K, T4coef(iF,:), T4corr(iF,:), T4pre(iF,2));
    fprintf(['                  | max rel. err coef: %.1e %.1e %.1e %.1e | corr: %.1e %.1e %.1e %.1e' ...
             ' | median angle err [rad]: grid %.1e, Lie-AD Newton %.1e\n'], ...
        E4coef(iF,:), E4corr(iF,:), median(E4ang(iF,:,4)), median(E4ang(iF,:,5)));
end

fig = new_fig(18, 8);
subplot(1,2,1);
for j = 1:4
    loglog(cfg.fixAxisL, T4coef(:,j)/K, ['-' mk{j}], 'Color', c(j,:), 'DisplayName', fa_methods{j}); hold on
end
set(gca, 'XScale', 'log', 'YScale', 'log'); grid on
xlabel('max band L'); ylabel('time per angle [s]'); legend('show', 'Location', 'northwest');
title(sprintf('Rotated coefficients, %d angles, fixed axis', K));
subplot(1,2,2);
for j = 1:4
    loglog(cfg.fixAxisL, T4corr(:,j)/K, ['-' mk{j}], 'Color', c(j,:), 'DisplayName', fa_methods{j}); hold on
end
set(gca, 'XScale', 'log', 'YScale', 'log'); grid on
xlabel('max band L'); ylabel('time per angle [s]'); legend('show', 'Location', 'northwest');
title(sprintf('Correlation C(\\theta), %d angles, fixed axis', K));
save_fig(fig, cfg, 'rsh_fixed_axis');


% ======================================================================
%% TEST 5: closed-form derivatives and Newton alignment on SO(3)
% ======================================================================
% For C(R) = sum_l d_l' M_l(R) s_l and the right perturbation R exp([w]x),
%   C = c + b'w + 1/2 w'Hw + O(|w|^3),  b_k = sum_l dt_l' L_k s_l,
%   H_jk = sum_l dt_l' (L_j L_k + L_k L_j) s_l / 2,  dt_l = M_l(R)' d_l,
% with the generators L_k (Section 5). (a) b and H are compared with central finite
% differences of C for several steps; (b) Newton iterations w = -H^-1 b, with a
% backtracking line search, recover the alignment of a vMF mixture from a 30-degree
% initial error.
L = cfg.gradL;
mus = rand_unit(numel(cfg.mix.w));
[axT, angT] = rand_axang(); Rtrue = axang_to_rotm(axT, angT);
s = mix_coeffs(L, mus, cfg.mix, eye(3), csConv);
d = mix_coeffs(L, mus, cfg.mix, Rtrue, csConv);
[ax0, ~] = rand_axang(); R0 = Rtrue * axang_to_rotm(ax0, 20*pi/180);
M0 = ad_all(L, rotm_to_rodrigues(R0), gens, sgnAD);
[b0, H0] = grad_hess(M0, s, d, gens, sgn, L);
C0 = corr_from_M(M0, s, d, L);
nH = numel(cfg.fdSteps); E5g = nan(1,nH); E5H = nan(1,nH);
for ih = 1:nH
    h = cfg.fdSteps(ih); E3 = eye(3); gfd = zeros(3,1); Hfd = zeros(3);
    Cw = @(w) corr_from_M(mult_M(M0, ad_all(L, w, gens, sgnAD), L), s, d, L);
    for k = 1:3
        Cp = Cw(h*E3(:,k)); Cm = Cw(-h*E3(:,k));
        gfd(k) = (Cp - Cm)/(2*h);
        Hfd(k,k) = (Cp - 2*C0 + Cm)/h^2;
        for j = k+1:3
            Hfd(j,k) = (Cw(h*(E3(:,j)+E3(:,k))) - Cw(h*(E3(:,j)-E3(:,k))) ...
                      - Cw(h*(E3(:,k)-E3(:,j))) + Cw(-h*(E3(:,j)+E3(:,k)))) / (4*h^2);
            Hfd(k,j) = Hfd(j,k);
        end
    end
    E5g(ih) = norm(gfd - b0)/norm(b0);
    E5H(ih) = norm(Hfd - H0, 'fro')/norm(H0, 'fro');
end
t5an = timeit(@() grad_hess(ad_all(L, rotm_to_rodrigues(R0), gens, sgnAD), s, d, gens, sgn, L));
t5fd = timeit(@() fd_grad_hess(R0, s, d, gens, sgnAD, L, 1e-4));
[~, ihg] = min(E5g); [~, ihH] = min(E5H);
fprintf(['Test 5a, L=%d | best FD rel. error: gradient %.1e (h=%.0e), Hessian %.1e (h=%.0e) | ' ...
         'time: closed form %.2e s, FD (19 evaluations) %.2e s\n'], ...
    L, E5g(ihg), cfg.fdSteps(ihg), E5H(ihH), cfg.fdSteps(ihH), t5an, t5fd);

E5N = nan(cfg.newtonTrials, cfg.newtonIters+1); D5 = nan(cfg.newtonTrials, 7);
for t = 1:cfg.newtonTrials
    mus = rand_unit(numel(cfg.mix.w));
    [axT, angT] = rand_axang(); Rtrue = axang_to_rotm(axT, angT);
    s = mix_coeffs(L, mus, cfg.mix, eye(3), csConv);
    d = mix_coeffs(L, mus, cfg.mix, Rtrue, csConv);
    [ax0, ~] = rand_axang(); R = Rtrue * axang_to_rotm(ax0, cfg.newtonInit);
    Mc = ad_all(L, rotm_to_rodrigues(R), gens, sgnAD);
    E5N(t,1) = rot_dist(R, Rtrue);
    for it = 1:cfg.newtonIters
        [b, H] = grad_hess(Mc, s, d, gens, sgn, L);
        % Newton direction; eigenvalues of H made negative, so that w is an ascent direction
        [Vh, Dh] = eig((H + H')/2);
        w = Vh * ((Vh' * b) ./ max(abs(diag(Dh)), 1e-12*max(abs(diag(Dh)))));
        if norm(w) > 0.2, w = 0.2 * w / norm(w); end
        % backtracking (Armijo) line search on C; the full step is accepted near the optimum
        C0 = corr_from_M(Mc, s, d, L);
        for ls = 1:30
            Mn = mult_M(Mc, ad_all(L, w, gens, sgnAD), L);
            if corr_from_M(Mn, s, d, L) >= C0 + 1e-4 * (b' * w), break; end
            w = w / 2;
        end
        Mc = Mn;
        if norm(w) > 0, R = R * axang_to_rotm(w/norm(w), norm(w)); end
        E5N(t,it+1) = rot_dist(R, Rtrue);
    end
    % Diagnostics at the last iterate: correlation reached relative to the optimum
    % C(R*) = sum_l ||d_l||^2, gradient norm, Hessian eigenvalues, and recent progress.
    [b, H] = grad_hess(Mc, s, d, gens, sgn, L);
    Copt = 0; for l = 0:L, Copt = Copt + d{l+1}' * d{l+1}; end
    D5(t,:) = [E5N(t,end), abs(E5N(t,end) - E5N(t,end-5)), corr_from_M(Mc, s, d, L) / Copt, ...
               norm(b) / Copt, sort(eig((H + H')/2))' / Copt];                 %#ok<SAGROW>
end
fprintf('Test 5b, L=%d, %d trials, initial error %.0f deg | median rotation error per iteration [rad]:\n   ', ...
    L, cfg.newtonTrials, cfg.newtonInit*180/pi);
fprintf(' %.1e', median(E5N, 1)); fprintf('\n');
fprintf('   converged (error < 1e-10) in %d of %d trials\n', sum(E5N(:,end) < 1e-10), cfg.newtonTrials);
% Classification of every trial. Quantities scaled by C(R*):
%   converged      : rotation error < 1e-10
%   local maximum  : error not changing, gradient ~ 0, all Hessian eigenvalues < 0, C < C(R*)
%   not converged  : otherwise (still progressing, or stalled away from a stationary point)
fprintf('Test 5b diagnostics (per trial, values relative to C(R*)):\n');
fprintf('  trial | final err [rad] | change last 5 it | C/C(R*)  | |grad|   | eig(H) (sorted)            | status\n');
status5 = cell(cfg.newtonTrials, 1);
for t = 1:cfg.newtonTrials
    q = D5(t,:);
    if q(1) < 1e-10
        status5{t} = 'converged';
    elseif q(2) < 1e-8 && q(4) < 1e-8 && all(q(5:7) < 0) && q(3) < 1 - 1e-8
        status5{t} = 'local maximum';
    else
        status5{t} = 'not converged';
    end
    fprintf('  %5d | %15.2e | %16.2e | %8.5f | %8.1e | %8.1e %8.1e %8.1e | %s\n', t, q, status5{t});
end
fprintf('   %d converged, %d at a local maximum, %d not converged\n', ...
    sum(strcmp(status5,'converged')), sum(strcmp(status5,'local maximum')), sum(strcmp(status5,'not converged')));

fig = new_fig(18, 8);
subplot(1,2,1);
loglog(cfg.fdSteps, E5g, '-o', 'DisplayName', 'gradient (FD vs closed form)'); hold on
loglog(cfg.fdSteps, E5H, '-s', 'DisplayName', 'Hessian (FD vs closed form)');
set(gca, 'XScale', 'log', 'YScale', 'log'); grid on
xlabel('finite-difference step h'); ylabel('relative difference'); legend('show', 'Location', 'south');
title('Closed-form derivatives vs finite differences');
subplot(1,2,2);
semilogy(0:cfg.newtonIters, E5N' + 1e-17, 'Color', [0.7 0.7 0.7]); hold on
semilogy(0:cfg.newtonIters, median(E5N,1) + 1e-17, 'k-o', 'LineWidth', 1.5);
set(gca, 'YScale', 'log'); grid on
xlabel('Newton iteration'); ylabel('rotation error [rad]');
title(sprintf('Newton alignment on SO(3), L=%d', L));
save_fig(fig, cfg, 'rsh_derivatives');


%% Summary and saving
fprintf('\n   l |  t_Ivanic  t_LieExp  t_LieEig   t_LieAD (all bands) | E_exp     E_eig     E_AD     (median)\n');
for il = 1:nL
    fprintf('%4d | %9.2e %9.2e %9.2e %9.2e             | %9.2e %9.2e %9.2e\n', cfg.lList(il), ...
        T_all(il,:), median(E(il,:,1)), median(E(il,:,2)), median(E(il,:,3)));
end
ensure_dir(cfg.resDir);
save(fullfile(cfg.resDir, 'rsh_rot_experiments_results.mat'), 'cfg', 'methods', 'T_all', 'T_band', 'E', 'E_expAD', 'O', ...
     'errA', 'errR', 'dnrm', 'sgn', 'sgnAD', 'ivMode', 'csConv', 'resExp', 'resAD', 'resIv', 'errF', 'errFo', ...
     'T4coef', 'T4corr', 'T4pre', 'E4coef', 'E4corr', 'E4ang', 'E5g', 'E5H', 't5an', 't5fd', 'E5N', 'D5', 'status5');

% ======================================================================
%% LOCAL FUNCTIONS
% ======================================================================

% ---- Timed loops -----------------------------------------------------------
function M = run_ivanic(l, Rl, mode)
    for k = 1:numel(Rl)
        M = ivanic_all(l, Rl{k}, mode);
    end
end

% Lie-exp / Lie-eig over the bands in gs (all bands 1..l, or band l only).
function M = run_lie(gs, Rl, Xl, method, sgn, conv)
    for k = 1:numel(Rl)
        if conv, x = rotm_to_rodrigues(Rl{k}); else, x = Xl(:,k); end
        for b = 1:numel(gs)
            A = lie_alg(gs{b}, sgn*x);
            if strcmp(method, 'exp'), M = expm(A); else, M = lie_eig(A); end
        end
    end
end

% Lie-AD: all bands through rsh_rot_lie_eig, or band l only through ad_band.
function M = run_ad(l, Rl, Xl, sgnAD, conv, allBands, gens)
    for k = 1:numel(Rl)
        if conv, x = rotm_to_rodrigues(Rl{k}); else, x = Xl(:,k); end
        if allBands, M = rsh_rot_lie_eig(l, sgnAD*x, gens); else, M = ad_band(l, sgnAD*x, gens{l}); end
    end
end

% All four matrices of band l for one rotation (order as in 'methods').
function Ms = all_methods(l, R, x, gens, sgn, sgnAD, ivMode)
    A  = lie_alg(gens{l}, sgn*x);
    Ms = {ivanic_band(l, R, ivMode), expm(A), lie_eig(A), rsh_rot_lie_eig(l, sgnAD*x, gens)};
end

% ---- Lie-algebraic methods ---------------------------------------------------
function A = lie_alg(g, x)
    A = x(1)*g.Lx + x(2)*g.Ly + x(3)*g.Lz;
end

% Exponential through the eigendecomposition. A is real skew-symmetric with distinct
% eigenvalues i*m*theta for theta ~= 0; the real part removes the imaginary round-off.
function M = lie_eig(A)
    [V, Dg] = eig(A);
    M = real((V * diag(exp(diag(Dg)))) / V);
end

% Anti-diagonal factorization of band l only: same algorithm as rsh_rot_lie_eig,
% without the lower bands.
function M = ad_band(l, x, G)
    theta = norm(x);
    if theta < 1e-14, M = eye(2*l+1); return; end
    Mz = eye(2*l+1);
    for m = 1:l
        i1 = l+1-m; i2 = l+1+m; cm = cos(m*theta); sm = sin(m*theta);
        Mz(i1,i1) = cm; Mz(i1,i2) = sm; Mz(i2,i1) = -sm; Mz(i2,i2) = cm;
    end
    F = rsh_rot_to_z(l, x/theta, G);
    M = F * Mz * F';
end

% ---- Ivanic-Ruedenberg recursion ---------------------------------------------
% Ivanic & Ruedenberg, J. Phys. Chem. 100, 6342 (1996); corrections in J. Phys. Chem. A
% 102, 9099 (1998). Vectorized over the column index n (no function handles in the
% inner loop, so that the timing reflects the algorithm and not the call overhead).
% mode: 1 M(R), 2 M(R'), 3 D M(R) D, 4 D M(R') D with D = diag((-1)^m).
function Ml = ivanic_band(l, R, mode)
    Ms = ivanic_all(l, R, mode);
    Ml = Ms{l};
end

function Ms = ivanic_all(l, R, mode)
    if nargin < 3, mode = 1; end
    if mode == 2 || mode == 4, R = R'; end
    Pm = [0 1 0; 0 0 1; 1 0 0];             % (x,y,z) -> (y,z,x): l = 1 real SH order
    R1 = Pm * R * Pm';
    Ms = cell(1, l);
    Ms{1} = R1;
    for ll = 2:l
        Ms{ll} = ir_next(ll, R1, Ms{ll-1});
    end
    if mode >= 3
        for ll = 1:l
            d = (-1).^(-ll:ll)';
            Ms{ll} = (d*d') .* Ms{ll};
        end
    end
end

function Rl = ir_next(l, R1, Rp)
    % Band l from band 1 (R1, indices -1..1) and band l-1 (Rp, indices -(l-1)..l-1).
    N = 2*l+1;
    % P{i+2}(a+l, b+l+1) = P(i, l, a, b), a = -(l-1)..l-1, b = -l..l
    P = cell(1,3);
    for i = -1:1
        ri = R1(i+2,:);                       % [r(i,-1) r(i,0) r(i,1)]
        Pi = zeros(2*l-1, N);
        Pi(:,2:N-1) = ri(2) * Rp;
        Pi(:,N) = ri(3)*Rp(:,end) - ri(1)*Rp(:,1);
        Pi(:,1) = ri(3)*Rp(:,1)   + ri(1)*Rp(:,end);
        P{i+2} = Pi;
    end
    n = -l:l;
    den = (l+n).*(l-n); den(abs(n) == l) = 2*l*(2*l-1);
    Rl = zeros(N);
    for m = -l:l
        am = abs(m); d = double(m == 0);
        row = zeros(1, N);
        if am < l                              % U term
            row = row + sqrt((l+m)*(l-m)./den) .* P{2}(m+l,:);
        end
        if m == 0                              % V term
            V = P{3}(1+l,:) + P{1}(-1+l,:);
        elseif m > 0
            V = P{3}(m-1+l,:)*sqrt(1+(m==1)) - P{1}(-m+1+l,:)*(1-(m==1));
        else
            V = P{3}(m+1+l,:)*(1-(m==-1)) + P{1}(-m-1+l,:)*sqrt(1+(m==-1));
        end
        row = row + 0.5*sqrt((1+d)*(l+am-1)*(l+am)./den)*(1-2*d) .* V;
        if am < l-1 && m ~= 0                  % W term (zero coefficient otherwise)
            if m > 0
                W = P{3}(m+1+l,:) + P{1}(-m-1+l,:);
            else
                W = P{3}(m-1+l,:) - P{1}(-m+1+l,:);
            end
            row = row - 0.5*sqrt((l-am-1)*(l-am)./den) .* W;
        end
        Rl(m+l+1,:) = row;
    end
end

% ---- Additional experiments (4-5) ---------------------------------------------
% vMF mixture coefficients of bands 0..L for the mean directions R*mus (cell, {l+1}).
function c = mix_coeffs(L, mus, mix, R, cs)
    c = cell(1, L+1); P = R * mus;
    for l = 0:L
        c{l+1} = rsh_eval_ref(l, P, cs) * (mix.w(:)' .* vmf_g(l, mix.kappa(:)'))';
    end
end

% Closed-form coefficients for all angles of a sweep about ax (cell of n_l x K).
function c = mix_coeffs_sweep(L, mus, mix, ax, th, cs)
    nK = size(mus,2); K = numel(th); P = zeros(3, nK*K);
    for k = 1:K, P(:, (k-1)*nK+(1:nK)) = axang_to_rotm(ax, th(k)) * mus; end
    c = cell(1, L+1);
    for l = 0:L
        n = 2*l+1; wg = mix.w(:)' .* vmf_g(l, mix.kappa(:)');
        Y = reshape(rsh_eval_ref(l, P, cs), n, nK, K);
        c{l+1} = reshape(sum(Y .* repmat(wg, [n 1 K]), 2), n, K);
    end
end

% Rotated coefficients for all angles about a fixed axis.
function out = fa_coef(method, L, s, ax, th, gens, sgn, sgnAD, ivMode)
    K = numel(th); out = cell(1, L+1); out{1} = repmat(s{1}, 1, K);
    switch method
        case 'Ivanic'
            for l = 1:L, out{l+1} = zeros(2*l+1, K); end
            for k = 1:K
                Ms = ivanic_all(L, axang_to_rotm(ax, th(k)), ivMode);
                for l = 1:L, out{l+1}(:,k) = Ms{l} * s{l+1}; end
            end
        case 'Lie-exp'
            for l = 1:L
                A = sgn * lie_alg(gens{l}, ax); o = zeros(2*l+1, K);
                for k = 1:K, o(:,k) = expm(th(k)*A) * s{l+1}; end
                out{l+1} = o;
            end
        case 'Lie-eig (reuse)'
            for l = 1:L
                [V, D] = eig(sgn * lie_alg(gens{l}, ax)); lam = diag(D);
                w = V \ s{l+1};
                out{l+1} = real(V * (exp(lam*th) .* repmat(w, 1, K)));
            end
        case 'Lie-AD (reuse)'
            for l = 1:L
                F = rsh_rot_to_z(l, sgnAD*ax, gens{l});
                out{l+1} = F * mz_apply_all(l, F' * s{l+1}, th);
            end
    end
end

% Correlation C(theta_k) = sum_l d_l' M_l(theta_k) s_l for all angles.
function C = fa_corr(method, L, s, d, ax, th, gens, sgn, sgnAD, ivMode)
    K = numel(th); C = d{1}' * s{1} * ones(1, K);
    switch method
        case 'Ivanic'
            for k = 1:K
                Ms = ivanic_all(L, axang_to_rotm(ax, th(k)), ivMode);
                for l = 1:L, C(k) = C(k) + d{l+1}' * (Ms{l} * s{l+1}); end
            end
        case 'Lie-exp'
            for l = 1:L
                A = sgn * lie_alg(gens{l}, ax);
                for k = 1:K, C(k) = C(k) + d{l+1}' * (expm(th(k)*A) * s{l+1}); end
            end
        case 'Lie-eig (reuse)'
            for l = 1:L
                [V, D] = eig(sgn * lie_alg(gens{l}, ax)); lam = diag(D);
                pw = (V.' * d{l+1}) .* (V \ s{l+1});
                C = C + real(sum(repmat(pw, 1, K) .* exp(lam*th), 1));
            end
        case 'Lie-AD (reuse)'
            [a0, A, B] = ad_trig_coeffs(L, s, d, sgnAD*ax, gens);
            C = trig_eval(a0, A, B, th);
    end
end

% M_lz(theta_k) q for all angles (columns).
function Q = mz_apply_all(l, q, th)
    Q = repmat(q, 1, numel(th));
    for m = 1:l
        i1 = l+1-m; i2 = l+1+m; cm = cos(m*th); sm = sin(m*th);
        Q(i1,:) =  cm*q(i1) + sm*q(i2);
        Q(i2,:) = -sm*q(i1) + cm*q(i2);
    end
end

% Coefficients of C(theta) = a0 + sum_m A_m cos(m theta) + B_m sin(m theta).
function [a0, A, B] = ad_trig_coeffs(L, s, d, axs, gens)
    a0 = d{1}' * s{1}; A = zeros(L,1); B = zeros(L,1);
    for l = 1:L
        F = rsh_rot_to_z(l, axs, gens{l});
        p = F' * d{l+1}; q = F' * s{l+1};
        m = (1:l)'; im = l+1-m; ip = l+1+m;
        a0 = a0 + p(l+1)*q(l+1);
        A(1:l) = A(1:l) + p(im).*q(im) + p(ip).*q(ip);
        B(1:l) = B(1:l) + p(im).*q(ip) - p(ip).*q(im);
    end
end

function C = trig_eval(a0, A, B, th)
    m = (1:numel(A))';
    C = a0 + A' * cos(m*th) + B' * sin(m*th);
end

% Newton iterations on the trigonometric polynomial, from th0.
function th = trig_argmax(A, B, th0)
    m = (1:numel(A))'; th = th0;
    for it = 1:10
        f1 = sum(m .* (-A.*sin(m*th) + B.*cos(m*th)));
        f2 = -sum(m.^2 .* (A.*cos(m*th) + B.*sin(m*th)));
        if f2 >= 0, break; end
        th = th - f1/f2;
    end
end

% Grid maximum refined by a parabola through the three neighbouring samples.
function th = grid_argmax(C, thK)
    K = numel(C); [~, k] = max(C); h = thK(2) - thK(1);
    ym = C(mod(k-2, K)+1); y0 = C(k); yp = C(mod(k, K)+1);
    den = ym - 2*y0 + yp; off = 0;
    if den < 0, off = 0.5*(ym - yp)/den; end
    th = thK(k) + off*h;
end

function e = ang_dist(a, b)
    e = abs(mod(a - b + pi, 2*pi) - pi);
end

% Rotation angle between two rotations, accurate also for small angles.
function e = rot_dist(R1, R2)
    e = 2*asin(min(1, norm(R1 - R2, 'fro') / (2*sqrt(2))));
end

% Per-axis precomputations (timed separately in Test 4).
function F = F_all_bands(L, axs, gens)
    F = cell(1, L);
    for l = 1:L, F{l} = rsh_rot_to_z(l, axs, gens{l}); end
end
function E = eig_all_bands(L, ax, gens, sgn)
    E = cell(2, L);
    for l = 1:L, [E{1,l}, E{2,l}] = eig(sgn * lie_alg(gens{l}, ax)); end
end

% All bands 0..L by Lie-AD (cell, {l+1}).
function M = ad_all(L, x, gens, sgnAD)
    [~, M] = rsh_rot_lie_eig(L, sgnAD*x, gens);
end
function M = mult_M(M1, M2, L)
    M = M1; for l = 1:L, M{l+1} = M1{l+1} * M2{l+1}; end
end
function C = corr_from_M(M, s, d, L)
    C = 0; for l = 0:L, C = C + d{l+1}' * (M{l+1} * s{l+1}); end
end

% Closed-form gradient and Hessian of C(R exp([w]x)) at w = 0.
function [b, H] = grad_hess(M, s, d, gens, sgn, L)
    b = zeros(3,1); H = zeros(3);
    for l = 1:L
        g = gens{l}; Lk = {sgn*g.Lx, sgn*g.Ly, sgn*g.Lz};
        dt = M{l+1}' * d{l+1};
        v = cell(1,3); u = cell(1,3);
        for k = 1:3, v{k} = Lk{k} * s{l+1}; u{k} = Lk{k}' * dt; end   % u_k' = dt' L_k
        for k = 1:3
            b(k) = b(k) + dt' * v{k};
            for j = 1:3
                H(j,k) = H(j,k) + (u{j}' * v{k} + u{k}' * v{j}) / 2;
            end
        end
    end
end

% Gradient and Hessian by central finite differences (19 evaluations of C).
function [g, H] = fd_grad_hess(R, s, d, gens, sgnAD, L, h)
    M0 = ad_all(L, rotm_to_rodrigues(R), gens, sgnAD);
    Cw = @(w) corr_from_M(mult_M(M0, ad_all(L, w, gens, sgnAD), L), s, d, L);
    E3 = eye(3); g = zeros(3,1); H = zeros(3); C0 = Cw([0;0;0]);
    for k = 1:3
        Cp = Cw(h*E3(:,k)); Cm = Cw(-h*E3(:,k));
        g(k) = (Cp - Cm)/(2*h); H(k,k) = (Cp - 2*C0 + Cm)/h^2;
        for j = k+1:3
            H(j,k) = (Cw(h*(E3(:,j)+E3(:,k))) - Cw(h*(E3(:,j)-E3(:,k))) ...
                    - Cw(h*(E3(:,k)-E3(:,j))) + Cw(-h*(E3(:,j)+E3(:,k)))) / (4*h^2);
            H(k,j) = H(j,k);
        end
    end
end

% ---- Figures: white background and PDF export -------------------------------
% Figure with white background and black axes, independent of the MATLAB theme
% (light or dark desktop). Size in centimetres.
function fig = new_fig(w, h)
    fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 w h], ...
                 'PaperUnits', 'centimeters', 'PaperSize', [w h], 'PaperPositionMode', 'auto');
    try, theme(fig, 'light'); catch, end            % R2025a+: force the light theme
    try
        set(fig, 'DefaultAxesColor', 'w', 'DefaultAxesXColor', 'k', 'DefaultAxesYColor', 'k', ...
            'DefaultAxesZColor', 'k', 'DefaultAxesGridColor', [0.15 0.15 0.15], ...
            'DefaultAxesFontSize', 9, 'DefaultTextColor', 'k', ...
            'DefaultLegendColor', 'w', 'DefaultLegendTextColor', 'k', 'DefaultLegendEdgeColor', [0.3 0.3 0.3]);
    catch
    end
end

% Create folder d (and its parents) if it does not exist. exist(d, 'dir')
% alone is not enough: for a relative name it also finds folders with that
% name on the MATLAB path, so the check is done on the absolute path.
function ensure_dir(d)
    if ~(numel(d) > 0 && (d(1) == '/' || d(1) == '\' || any(d == ':')))
        d = fullfile(pwd, d);                          % make it absolute
    end
    if exist(d, 'dir') ~= 7
        [ok, msg] = mkdir(d);
        if ~ok
            error('Cannot create output folder %s: %s', d, msg);
        end
        fprintf('Created output folder %s\n', d);
    end
end

% Applies the white style to all axes and legends of fig and saves it as
% <cfg.figDir>/<name>.pdf (vector graphics, white background).
function save_fig(fig, cfg, name)
    ax = findall(fig, 'Type', 'axes');
    set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.15 0.15 0.15], 'Box', 'on');
    for a = ax'
        set(get(a, 'Title'), 'Color', 'k'); set(get(a, 'XLabel'), 'Color', 'k'); set(get(a, 'YLabel'), 'Color', 'k');
    end
    set(findall(fig, 'Type', 'legend'), 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.3 0.3 0.3]);
    if ~cfg.saveFigs, return; end
    ensure_dir(cfg.figDir);
    file = fullfile(cfg.figDir, [name '.pdf']);
    drawnow;
    if exist('exportgraphics', 'file')                 % R2020a+
        exportgraphics(fig, file, 'ContentType', 'vector', 'BackgroundColor', 'white');
    else
        set(fig, 'InvertHardcopy', 'off');
        print(fig, file, '-dpdf', '-painters');
    end
    fprintf('Saved %s\n', file);
end

% ---- Real spherical harmonics and vMF ----------------------------------------
% Orthonormal real SH of band l at the unit vectors U (3 x N), rows m = -l..l,
% m < 0: sin(|m| phi), m > 0: cos(m phi). cs = 1 adds the phase (-1)^m.
% legendre(...,'norm') does not include the Condon-Shortley phase.
function Y = rsh_eval_ref(l, U, cs)
    z = max(-1, min(1, U(3,:))); phi = atan2(U(2,:), U(1,:));
    P = legendre(l, z, 'norm');               % (l+1) x N, rows m = 0..l
    if l == 0, P = reshape(P, 1, []); end
    Y = zeros(2*l+1, size(U,2));
    Y(l+1,:) = P(1,:) / sqrt(2*pi);
    for m = 1:l
        cm = P(m+1,:) / sqrt(pi);
        if cs, cm = (-1)^m * cm; end
        Y(l+1+m,:) = cm .* cos(m*phi);
        Y(l+1-m,:) = cm .* sin(m*phi);
    end
end

function r = equiv_res(M, l, R, U, cs)
    r = max(sqrt(sum((rsh_eval_ref(l, R*U, cs) - M*rsh_eval_ref(l, U, cs)).^2, 1)));
end

% g_l(kappa) = kappa i_l(kappa) / sinh(kappa), i_l(k) = sqrt(pi/(2k)) I_{l+1/2}(k),
% evaluated with the exponentially scaled Bessel function. l and kappa may be vectors
% (l column-wise broadcast against kappa row-wise when both are vectors).
function g = vmf_g(l, kappa)
    [L, K] = ndgrid(l(:), kappa(:));
    g = sqrt(pi*K/2) .* 2 .* besseli(L + 0.5, K, 1) ./ (1 - exp(-2*K));
    if isscalar(l), g = reshape(g, 1, []); elseif isscalar(kappa), g = reshape(g, 1, []); end
end

% ---- Rotations ---------------------------------------------------------------
% Haar-distributed random rotation as axis-angle (angle in [0, pi]).
function [ax, ang] = rand_axang()
    q = randn(4,1); q = q/norm(q); if q(1) < 0, q = -q; end
    ang = 2*acos(min(1, q(1)));
    v = q(2:4);
    if norm(v) < eps, ax = [0;0;1]; else, ax = v/norm(v); end
end

function U = rand_unit(n)
    U = randn(3, n); U = U ./ sqrt(sum(U.^2, 1));
end

function R = axang_to_rotm(ax, ang)
    K = [0 -ax(3) ax(2); ax(3) 0 -ax(1); -ax(2) ax(1) 0];
    R = eye(3) + sin(ang)*K + (1-cos(ang))*(K*K);
end

% Rodrigues vector from a rotation matrix, robust near 0 and pi.
function x = rotm_to_rodrigues(R)
    c = max(-1, min(1, (trace(R) - 1)/2));
    ang = acos(c);
    v = [R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)];
    if ang < 1e-12
        x = v/2;                                 % first-order, exact limit at 0
    elseif pi - ang > 1e-6
        x = ang * v / (2*sin(ang));
    else                                         % near pi: axis from R + I
        B = (R + eye(3))/2;
        [~, i] = max(diag(B));
        ax = B(:,i) / sqrt(B(i,i));
        if ax' * v < 0, ax = -ax; end
        x = ang * ax / norm(ax);
    end
end
