function run_tests()
%RUN_TESTS  Self-checks of the library (MATLAB R2016b+ or GNU Octave 6+).
%
%   Run from any folder:
%       run_tests
%   Every check prints PASS/FAIL with the measured error; the function
%   throws an error at the end if any check failed. Takes a few seconds.

  repoDir = fileparts(fileparts(mfilename('fullpath')));
  if exist(fullfile(repoDir, 'src', 'rsh_lie_generator.m'), 'file') ~= 2
    repoDir = fileparts(pwd);                    % fallback: run from tests/
  end
  addpath(genpath(fullfile(repoDir, 'src')));
  rng(1);
  nFail = 0;
  tol = 1e-10;

  % ---- 1. Lie algebra generators -------------------------------------------
  e = 0;
  for l = 1:12
    [Lx, Ly, Lz] = rsh_lie_generator(l);
    I = eye(2*l+1);
    e = max([e, norm(Lx + Lx', 'fro'), norm(Ly + Ly', 'fro'), norm(Lz + Lz', 'fro'), ...
             norm(Lx*Ly - Ly*Lx - Lz, 'fro'), norm(Ly*Lz - Lz*Ly - Lx, 'fro'), ...
             norm(Lz*Lx - Lx*Lz - Ly, 'fro'), ...
             norm(Lx^2 + Ly^2 + Lz^2 + l*(l+1)*I, 'fro') / (l*(l+1))]);
  end
  nFail = nFail + check('generators: skew-symmetry, [Lx,Ly]=Lz, Casimir', e, tol);

  % ---- 2. Anti-diagonal factorization vs matrix exponential -----------------
  lmax = 20;
  X = [randn(3, 20), [0; 0; 1.3], [0; 0; -0.4], 1e-13 * [1; 1; 1], ...
       pi * [1; 0; 0], (pi - 1e-9) * [0; 1; 1] / sqrt(2)];
  nFail = nFail + check('rsh_rot_lie_eig == expm (incl. z axis, tiny and pi angles)', ...
                        max_diff_ad_exp(lmax, X), 1e-11);
  % Known limitation (see the note in rsh_rot_to_z): for an axis at angle rho
  % from +/-z (0 < rho << 1) the accuracy of rsh_eval limits F_l to ~eps/rho.
  nFail = nFail + check('rsh_rot_lie_eig, axis 2e-9 rad from z (error ~eps/rho)', ...
                        max_diff_ad_exp(lmax, [1e-9; 2e-9; 0.7]), 1e-7);

  % ---- 3. Orthogonality and F_l self-check ------------------------------------
  e = 0; eF = 0;
  for k = 1:5
    x = randn(3, 1);
    [~, M] = rsh_rot_lie_eig(32, x);
    for l = [1 8 16 32]
      e = max(e, norm(M{l+1}' * M{l+1} - eye(2*l+1), 'fro'));
      [Lx, Ly, Lz] = rsh_lie_generator(l);
      F = rsh_rot_to_z(l, x);
      u = x / norm(x);
      eF = max([eF, norm(F' * F - eye(2*l+1), 'fro'), ...
                norm(F' * (u(1)*Lx + u(2)*Ly + u(3)*Lz) * F - Lz, 'fro') / l]);
    end
  end
  nFail = nFail + check('M_l orthogonal (l up to 32)', e, 1e-11);
  nFail = nFail + check('rsh_rot_to_z: F orthogonal, F''*(u.L)*F = Lz', eF, 1e-11);

  % ---- 4. Equivariance Y_l(R u) = M_l Y_l(u) -----------------------------------
  lmax = 10; e = 0;
  for k = 1:10
    x = randn(3, 1); R = exp_SO3(x);
    u = randn(3, 1); u = u / norm(u);
    Y = rsh_eval(lmax, u); Yr = rsh_eval(lmax, R * u);
    [~, M] = rsh_rot_lie_eig(lmax, x);
    for l = 0:lmax
      idx = sh_lm_to_index(l, -l:l);
      e = max(e, norm(Yr(idx) - M{l+1} * Y(idx)));
    end
  end
  nFail = nFail + check('equivariance Y(R u) = M(R) Y(u)', e, tol);

  % ---- 5. Group homomorphism M(R1 R2) = M(R1) M(R2) ----------------------------
  e = 0;
  for k = 1:10
    x1 = randn(3, 1); x2 = randn(3, 1);
    x12 = log_SO3(exp_SO3(x1) * exp_SO3(x2));
    [~, M1] = rsh_rot_lie_eig(8, x1); [~, M2] = rsh_rot_lie_eig(8, x2);
    [~, M12] = rsh_rot_lie_eig(8, x12);
    for l = 1:8, e = max(e, norm(M12{l+1} - M1{l+1} * M2{l+1}, 'fro')); end
  end
  nFail = nFail + check('homomorphism M(R1 R2) = M(R1) M(R2)', e, 1e-9);

  % ---- 6. Ivanic-Ruedenberg baseline ------------------------------------------
  e = 0;
  for k = 1:5
    x = randn(3, 1);
    [~, Mi] = rsh_rot_ivanic(exp_SO3(x), 6);
    [~, Ml] = rsh_rot_lie_eig(6, log_SO3(exp_SO3(x)));
    for l = 0:6, e = max(e, norm(Mi{l+1} - Ml{l+1}, 'fro')); end
  end
  nFail = nFail + check('rsh_rot_ivanic(R) == rsh_rot_lie_eig(log R)', e, 1e-10);

  % ---- 7. Rotation of an expansion (closed-form vMF coefficients) ------------
  lmax = 24; kappa = 15; mu = [1; -2; 2] / 3; x = [0.4; 1.1; -0.6];
  c   = von_mises_fisher_rsh_coeffs(lmax, mu, kappa);
  cr  = rsh_rot_coeffs(x, c);
  cr2 = rsh_rot_coeffs(exp_SO3(x), c);
  ref = von_mises_fisher_rsh_coeffs(lmax, exp_SO3(x) * mu, kappa);
  nFail = nFail + check('rsh_rot_coeffs on vMF (Rodrigues input)', norm(cr - ref) / norm(ref), tol);
  nFail = nFail + check('rsh_rot_coeffs on vMF (matrix input)', norm(cr2 - ref) / norm(ref), tol);

  % ---- 8. Orthonormality of rsh_eval / csh_eval (midpoint quadrature) --------
  n = 200; lmax = 3;
  th = ((1:n) - 0.5) * pi / n; ph = ((1:2*n) - 0.5) * pi / n;
  w = sin(th(:)) * (pi / n)^2;
  Yr = reshape(rsh_eval(lmax, th, ph), (lmax+1)^2, []);
  Yc = reshape(csh_eval(lmax, th, ph), (lmax+1)^2, []);
  W = repmat(w, 2*n, 1).';
  Gr = (Yr .* W) * Yr'; Gc = (Yc .* W) * Yc';
  nFail = nFail + check('rsh_eval orthonormal', norm(Gr - eye(size(Gr)), 'fro'), 1e-3);
  nFail = nFail + check('csh_eval orthonormal', norm(Gc - eye(size(Gc)), 'fro'), 1e-3);

  % ---- 9. Correlation: closed-form gradient / Hessian vs finite differences --
  lmax = 6; s = randn((lmax+1)^2, 1); d = randn((lmax+1)^2, 1);
  R0 = exp_SO3([0.3; -0.7; 0.5]);
  [f0, g, H] = rsh_corr_grad_hess(lmax, s, d, R0);
  fun = @(w) rsh_corr_grad_hess(lmax, s, d, R0 * exp_SO3(w));
  h = 1e-4; E = eye(3); gfd = zeros(3, 1); Hfd = zeros(3);
  for j = 1:3
    gfd(j) = (fun(h*E(:,j)) - fun(-h*E(:,j))) / (2*h);
    for k = 1:3
      Hfd(j,k) = (fun(h*(E(:,j)+E(:,k))) - fun(h*(E(:,j)-E(:,k))) ...
                - fun(h*(E(:,k)-E(:,j))) + fun(-h*(E(:,j)+E(:,k)))) / (4*h^2);
    end
  end
  nFail = nFail + check('rsh_corr_grad_hess gradient vs FD', norm(g - gfd) / norm(g), 1e-6);
  nFail = nFail + check('rsh_corr_grad_hess Hessian vs FD', norm(H - Hfd, 'fro') / norm(H, 'fro'), 1e-5);

  % ---- 10. Newton alignment on SO(3) ------------------------------------------
  lmax = 8; mu = [0.2; 0.5; 0.84]; mu = mu / norm(mu);
  s = von_mises_fisher_rsh_coeffs(lmax, mu, 10) + 0.5 * von_mises_fisher_rsh_coeffs(lmax, [1; 0; 0], 30);
  Rtrue = exp_SO3([0.5; -0.3; 0.8]);
  d = rsh_rot_coeffs(Rtrue, s);
  Rinit = Rtrue * exp_SO3(0.3 * [1; 1; -1] / sqrt(3));     % ~17 degrees away
  Ropt = rsh_opt_newton_so3(lmax, s, d, Rinit, struct('verbose', false, 'maxIter', 50));
  nFail = nFail + check('rsh_opt_newton_so3 recovers the rotation', norm(log_SO3(Rtrue' * Ropt)), 1e-6);

  % ---- 11. SO(3) utilities ---------------------------------------------------------
  e = 0;
  for k = 1:20
    x = randn(3, 1); x = x / norm(x) * (pi * rand);
    e = max(e, norm(log_SO3(exp_SO3(x)) - x));
  end
  nFail = nFail + check('log_SO3(exp_SO3(x)) = x', e, 1e-9);
  a = [0.3; -0.2; 0.9]; b = [-0.5; 0.4; 0.1]; l = 4;
  [Om, A, B] = adjoint_rodrigues(a, b, l);
  nFail = nFail + check('adjoint_rodrigues = expm(A) B expm(-A)', ...
                        norm(Om - expm(A) * B * expm(-A), 'fro'), 1e-12);

  % ---- summary -------------------------------------------------------------------
  if nFail > 0
    error('run_tests:failed', '%d check(s) FAILED.', nFail);
  end
  fprintf('All checks passed.\n');
end

function e = max_diff_ad_exp(lmax, X)
  e = 0;
  for k = 1:size(X, 2)
    [~, Mad]  = rsh_rot_lie_eig(lmax, X(:,k));
    [~, Mexp] = rsh_rot_lie(lmax, X(:,k));
    for l = 0:lmax
      e = max(e, norm(Mad{l+1} - Mexp{l+1}, 'fro') / (2*l+1));
    end
  end
end

function failed = check(name, err, tol)
  failed = ~(err <= tol);
  if failed, status = 'FAIL'; else, status = 'PASS'; end
  fprintf('[%s] %-62s err = %.2e (tol %.0e)\n', status, name, err, tol);
end
