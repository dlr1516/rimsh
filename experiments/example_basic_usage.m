%% EXAMPLE_BASIC_USAGE  Minimal tour of the library (runs in a few seconds).

repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repoDir, 'src')));

%% 1. A function on the sphere: von Mises-Fisher density, bands 0..lmax
lmax  = 16;
mu    = [0; 0; 1];                         % mean direction
kappa = 30;                                % concentration
c = von_mises_fisher_rsh_coeffs(lmax, mu, kappa);   % (lmax+1)^2 RSH coefficients

%% 2. Rotation matrices of all bands with the proposed method
x = (pi/3) * [1; 1; 0] / sqrt(2);          % Rodrigues vector: 60 deg about [1 1 0]
[~, M] = rsh_rot_lie_eig(lmax, x);         % M{l+1} is the (2l+1)x(2l+1) matrix of band l

% Same result with the generic matrix exponential of the Lie algebra element
[Lx, Ly, Lz] = rsh_lie_generator(lmax);
Mexp = expm(x(1)*Lx + x(2)*Ly + x(3)*Lz);
fprintf('||M_lmax(anti-diagonal) - expm||_F = %.2e\n', norm(M{end} - Mexp, 'fro'));

%% 3. Rotate the expansion and check it against the rotated vMF
c_rot = rsh_rot_coeffs(x, c);
c_ref = von_mises_fisher_rsh_coeffs(lmax, exp_SO3(x) * mu, kappa);
fprintf('||c_rot - c_ref|| / ||c_ref|| = %.2e\n', norm(c_rot - c_ref) / norm(c_ref));

%% 4. Plot before / after
figure;
subplot(1,2,1); rsh_plot(c,     90, 'vMF, original');
subplot(1,2,2); rsh_plot(c_rot, 90, 'vMF, rotated');
