%% Timing and accuracy of RSH rotation methods vs band l (DRAFT, untested)
% Methods (all produce the (2l+1)x(2l+1) matrix M_l of ONE band l):
%   'ivanic'   Ivanic-Ruedenberg recursion (reference). It needs all bands 1..l, so
%              its cost for band l includes computing the lower bands.
%   'lie-exp'  M_l = expm(theta * n.L^l)          (generic matrix exponential)
%   'lie-eig'  M_l = V * exp(Lambda) * V^-1 with [V,Lambda] = eig(theta * n.L^l)
%
% Requires: rsh_lie_generator(l) (single band, size (2l+1)x(2l+1)), axang2rotm, rotm2axang.
% Test 1: computation time per rotation vs l.
% Test 2: difference of Lie-exp and Lie-eig from Ivanic vs l (and Lie-exp vs Lie-eig).

%% 0. Configuration
clear; clc; rng(3);
cfg.lList      = [1 2 3 4 6 8 12 16 24 32 48 64];   % bands to test
cfg.nRotTime   = 100;     % rotations per timing measurement (time reported per rotation)
cfg.nRotAcc    = 50;      % random rotations per l for the accuracy test
cfg.includeConversion = true;   % true: every method receives a 3x3 matrix R, so the Lie
                                % methods also pay for rotm2axang (fair comparison);
                                % false: Lie methods receive axis and angle directly
Pm = [0 1 0; 0 0 1; 1 0 0];     % (x,y,z) -> (y,z,x), order of the l = 1 real SH
nL = numel(cfg.lList);

% Generators are computed once per l and are NOT part of the timings
% (they do not depend on the rotation; mention this in the paper).
gens = cell(1, nL);
for il = 1:nL
    l = cfg.lList(il);
    [gx, gy, gz] = rsh_lie_generator(l);
    assert(size(gx,1) == 2*l+1, 'rsh_lie_generator(%d) has size %d', l, size(gx,1));
    gens{il} = struct('Lx', gx, 'Ly', gy, 'Lz', gz);
end

%% 1. Convention check between generators and the recursion
% Ivanic is the reference, so the generators must use the same basis and sign
% convention. For l = 2 both signs of the exponent are compared with Ivanic; the
% better one is used from here on (sgn = +1 means M = exp(+theta n.L)).
% A large error for both signs means the SH basis conventions differ.
g2 = struct('Lx', [], 'Ly', [], 'Lz', []);
[g2.Lx, g2.Ly, g2.Lz] = rsh_lie_generator(2);
ep = 0; em = 0;
for k = 1:20
    [ax, ang] = rand_axang(); R = axang2rotm([ax' ang]);
    Mi = ir_last_band(2, R);
    A  = ang*(ax(1)*g2.Lx + ax(2)*g2.Ly + ax(3)*g2.Lz);
    ep = max(ep, norm(expm( A) - Mi,'fro'));
    em = max(em, norm(expm(-A) - Mi,'fro'));
end
if ep <= em, sgn = +1; else, sgn = -1; end
fprintf('Convention check (l=2): error with +A = %.2e, with -A = %.2e -> sgn = %+d\n', ep, em, sgn);
if min(ep, em) > 1e-8
    warning(['Neither sign matches Ivanic: basis convention differs or the recursion ' ...
             'implementation is wrong. Results below are not meaningful until this is fixed.']);
end

% ======================================================================
%% TEST 1: computation time vs l
% ======================================================================
% For each l the SAME list of random rotations is used for all methods. timeit
% runs warm-up and repeats the whole loop of cfg.nRotTime rotations; the reported
% time is per rotation.
T = nan(nL, 3);                              % columns: ivanic, lie-exp, lie-eig
for il = 1:nL
    l = cfg.lList(il); g = gens{il};
    Rl = cell(cfg.nRotTime,1); AXl = zeros(cfg.nRotTime, 4);
    for k = 1:cfg.nRotTime
        [ax, ang] = rand_axang();
        Rl{k} = axang2rotm([ax' ang]); AXl(k,:) = [ax' ang];
    end
    T(il,1) = timeit(@() run_ivanic(l, Rl))                           / cfg.nRotTime;
    T(il,2) = timeit(@() run_lie(g, Rl, AXl, 'exp', sgn, cfg.includeConversion)) / cfg.nRotTime;
    T(il,3) = timeit(@() run_lie(g, Rl, AXl, 'eig', sgn, cfg.includeConversion)) / cfg.nRotTime;
    fprintf('l=%2d | time per rotation [s]: ivanic %.2e | lie-exp %.2e | lie-eig %.2e\n', l, T(il,:));
end
figure;
loglog(cfg.lList, T, '-o', 'LineWidth', 1.2); grid on
xlabel('band l'); ylabel('time per rotation [s]');
legend('Ivanic', 'Lie-exp', 'Lie-eig', 'Location', 'northwest');
title('Computation time of M_l vs l');

% ======================================================================
%% TEST 2: difference from Ivanic (reference) vs l
% ======================================================================
% For every l the same random rotations are processed by the three methods.
% Errors are Frobenius norms of the difference of the (2l+1)x(2l+1) matrices.
% E_ei = Lie-exp  vs Ivanic,   E_gi = Lie-eig vs Ivanic,   E_eg = Lie-exp vs Lie-eig
E_ei = nan(nL, cfg.nRotAcc); E_gi = E_ei; E_eg = E_ei;
for il = 1:nL
    l = cfg.lList(il); g = gens{il};
    for k = 1:cfg.nRotAcc
        [ax, ang] = rand_axang(); R = axang2rotm([ax' ang]);
        Mi = ir_last_band(l, R);
        A  = sgn*ang*(ax(1)*g.Lx + ax(2)*g.Ly + ax(3)*g.Lz);
        Me = expm(A);
        Mg = lie_eig(A);
        E_ei(il,k) = norm(Me - Mi,'fro');
        E_gi(il,k) = norm(Mg - Mi,'fro');
        E_eg(il,k) = norm(Me - Mg,'fro');
    end
    fprintf('l=%2d | median |M-Ivanic|_F: lie-exp %.2e | lie-eig %.2e | lie-exp vs lie-eig %.2e\n', ...
        l, median(E_ei(il,:)), median(E_gi(il,:)), median(E_eg(il,:)));
end
figure; hold on
c = lines(3);
semilogy(cfg.lList, median(E_ei,2) + eps, '-o',  'Color', c(1,:), 'DisplayName', 'Lie-exp vs Ivanic (median)');
semilogy(cfg.lList, max(E_ei,[],2)  + eps, '--o', 'Color', c(1,:), 'DisplayName', 'Lie-exp vs Ivanic (max)');
semilogy(cfg.lList, median(E_gi,2) + eps, '-s',  'Color', c(2,:), 'DisplayName', 'Lie-eig vs Ivanic (median)');
semilogy(cfg.lList, max(E_gi,[],2)  + eps, '--s', 'Color', c(2,:), 'DisplayName', 'Lie-eig vs Ivanic (max)');
semilogy(cfg.lList, median(E_eg,2) + eps, '-^',  'Color', c(3,:), 'DisplayName', 'Lie-exp vs Lie-eig (median)');
set(gca, 'YScale', 'log', 'XScale', 'log'); grid on
xlabel('band l'); ylabel('||M - M_{ref}||_F'); legend('show', 'Location', 'northwest');
title('Difference from Ivanic vs l');

%% Summary table and saving
Tab = table(cfg.lList(:), T(:,1), T(:,2), T(:,3), median(E_ei,2), max(E_ei,[],2), ...
    median(E_gi,2), max(E_gi,[],2), median(E_eg,2), ...
    'VariableNames', {'l','t_ivanic','t_lie_exp','t_lie_eig','err_exp_med','err_exp_max', ...
                      'err_eig_med','err_eig_max','diff_exp_eig_med'});
disp(Tab);
save('rsh_rot_timing_accuracy_results.mat', 'cfg', 'T', 'E_ei', 'E_gi', 'E_eg', 'Tab', 'sgn');

% ======================================================================
%% LOCAL FUNCTIONS
% ======================================================================

% ---- Loops used by timeit --------------------------------------------------
% Ivanic: input is the 3x3 matrix; all bands up to l are built, the last one is kept.
function M = run_ivanic(l, Rl)
    for k = 1:numel(Rl)
        M = ir_last_band(l, Rl{k});
    end
end

% Lie methods: 'exp' or 'eig'. With conv = true the axis-angle is recovered from R
% inside the timed loop, otherwise it is taken from the precomputed list AXl.
function M = run_lie(g, Rl, AXl, method, sgn, conv)
    for k = 1:numel(Rl)
        if conv
            a4 = rotm2axang(Rl{k});
        else
            a4 = AXl(k,:);
        end
        A = sgn * a4(4) * (a4(1)*g.Lx + a4(2)*g.Ly + a4(3)*g.Lz);
        if strcmp(method, 'exp')
            M = expm(A);
        else
            M = lie_eig(A);
        end
    end
end

% ---- Lie-eig: exponential through the eigendecomposition --------------------
% A is real skew-symmetric (normal) with eigenvalues i*m*theta, m = -l..l, which are
% distinct for theta ~= 0. The result is real up to round-off, so the real part is taken.
function M = lie_eig(A)
    [V, Dg] = eig(A);
    M = real((V * diag(exp(diag(Dg)))) / V);
end

% ---- Ivanic-Ruedenberg recursion ---------------------------------------------
% Last band (size (2l+1)x(2l+1)) of the recursion started from D^1 = P R P'.
% Reference: Ivanic & Ruedenberg, J. Phys. Chem. 100, 6342 (1996), with the
% corrections of J. Phys. Chem. A 102, 9099 (1998). Formulas typed from memory:
% if the convention check at the top reports a mismatch, verify them first.
function Ml = ir_last_band(l, R)
    Pm = [0 1 0; 0 0 1; 1 0 0];
    R1 = Pm * R * Pm';                      % band 1, indices m = -1, 0, 1
    if l == 0, Ml = 1; return; end
    Prev = R1;
    for ll = 2:l
        Prev = ir_next(ll, R1, Prev);
    end
    Ml = Prev;
end

function Rl = ir_next(l, R1, Rp)
    % Band l from band l-1 (Rp) and band 1 (R1); indices m, n in -l..l
    Rl = zeros(2*l+1);
    r  = @(i,j) R1(i+2, j+2);               % R^1 with indices -1..1
    rp = @(a,b) Rp(a+l, b+l);               % R^{l-1} with indices -(l-1)..(l-1)
    for m = -l:l
        for n = -l:l
            d = double(m == 0);
            if abs(n) < l, den = (l+n)*(l-n); else, den = 2*l*(2*l-1); end
            u = sqrt((l+m)*(l-m)/den);
            v = 0.5*sqrt((1+d)*(l+abs(m)-1)*(l+abs(m))/den)*(1-2*d);
            w = -0.5*sqrt((l-abs(m)-1)*(l-abs(m))/den)*(1-d);
            val = 0;
            % Terms with zero coefficient are skipped (their indices may be out of range)
            if u ~= 0, val = val + u*Pf(0, m, n); end
            if v ~= 0, val = val + v*Vf(m, n);    end
            if w ~= 0, val = val + w*Wf(m, n);    end
            Rl(m+l+1, n+l+1) = val;
        end
    end

    function val = Pf(i, mu, nn)
        if nn == l
            val = r(i,1)*rp(mu, l-1) - r(i,-1)*rp(mu, -l+1);
        elseif nn == -l
            val = r(i,1)*rp(mu, -l+1) + r(i,-1)*rp(mu, l-1);
        else
            val = r(i,0)*rp(mu, nn);
        end
    end
    function val = Vf(mm, nn)
        if mm == 0
            val = Pf(1,1,nn) + Pf(-1,-1,nn);
        elseif mm > 0
            val = Pf(1,mm-1,nn)*sqrt(1+(mm==1)) - Pf(-1,-mm+1,nn)*(1-(mm==1));
        else
            val = Pf(1,mm+1,nn)*(1-(mm==-1)) + Pf(-1,-mm-1,nn)*sqrt(1+(mm==-1));
        end
    end
    function val = Wf(mm, nn)
        if mm > 0
            val = Pf(1,mm+1,nn) + Pf(-1,-mm-1,nn);
        else
            val = Pf(1,mm-1,nn) - Pf(-1,-mm+1,nn);
        end
    end
end

% ---- Haar-distributed random rotation as axis-angle ---------------------------
function [ax, ang] = rand_axang()
    q = randn(4,1); q = q/norm(q); if q(1) < 0, q = -q; end
    ang = 2*acos(min(1, q(1)));
    v = q(2:4);
    if norm(v) < eps, ax = [0;0;1]; else, ax = v/norm(v); end
end