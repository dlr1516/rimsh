function demo_von_mises_fisher(lmax, mu, kappa, res)
%DEMO_VON_MISES_FISHER  Compare a von Mises-Fisher function with its RSH series.
%
%   demo_von_mises_fisher(lmax, mu, kappa, res)
%
%   Defaults: lmax = 5, mu = unit vector [1;2;2]/3, kappa = 1.0, res = pi/180.
%   Needs: rsh_eval.m, rsh_plot.m, von_mises_fisher_rsh_coeffs.m (src/ is added
%   to the path automatically).

repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repoDir, 'src')));

if nargin < 1 || isempty(lmax),  lmax  = 5; end
if nargin < 2 || isempty(mu),    mu    = [1; 2; 2] / 3; end
if nargin < 3 || isempty(kappa), kappa = 1.0; end
if nargin < 4 || isempty(res),   res   = pi/180; end

mu = mu(:) / norm(mu);

theta = 0:res:pi;            % colatitude
phi   = 0:res:(2*pi);        % azimuth
nt = length(theta);
np = length(phi);

% ---- 1) exact vMF function on the grid ----
[Phi, Theta] = meshgrid(phi, theta);                 % both nt x np
ux = sin(Theta).*cos(Phi);
uy = sin(Theta).*sin(Phi);
uz = cos(Theta);
dotmu = mu(1)*ux + mu(2)*uy + mu(3)*uz;

% f = kappa/(4*pi*sinh(kappa)) * exp(kappa*dotmu), written to avoid overflow
if kappa == 0
    F = ones(nt, np) / (4*pi);
else
    F = kappa / (2*pi*(-expm1(-2*kappa))) * exp(kappa*(dotmu - 1));
end

% ---- 2) RSH coefficients and series expansion ----
c = von_mises_fisher_rsh_coeffs(lmax, mu, kappa);    % (lmax+1)^2 x 1
Y = rsh_eval(lmax, theta, phi);                      % (lmax+1)^2 x nt x np
Fseries = reshape(c(:).' * reshape(Y, [], nt*np), nt, np);

err = Fseries - F;
fprintf('c(1) = %.10f   (1/sqrt(4*pi) = %.10f)\n', c(1), 1/sqrt(4*pi));
fprintf('max |f|        = %.6g\n', max(abs(F(:))));
fprintf('max |f - f_L|  = %.6g  (lmax = %d)\n', max(abs(err(:))), lmax);
fprintf('rel. L2 error  = %.6g\n', norm(err(:)) / norm(F(:)));

latDeg = 90 - 180/pi*theta;
lonDeg = 180/pi*phi;

% ---- 3) plots ----
figure;

subplot(2,2,1);
pc = pcolor(lonDeg, latDeg, F);
set(pc, 'EdgeColor', 'none');
title('von Mises-Fisher f(\theta,\phi)');
xlabel('lon [deg]'); ylabel('lat [deg]'); colorbar();

subplot(2,2,2);
% rsh_plot uses step = pi/res_plot, so res_plot = pi/res
rsh_plot(c, round(pi/res), sprintf('RSH series, lmax = %d', lmax), 'none');

subplot(2,2,3);
pc = pcolor(lonDeg, latDeg, err);
set(pc, 'EdgeColor', 'none');
title('Series - exact');
xlabel('lon [deg]'); ylabel('lat [deg]'); colorbar();

subplot(2,2,4);
stem(0:(lmax+1)^2-1, c, 'filled');
title('RSH coefficients (packed index l^2+l+m)');
xlabel('index - 1'); ylabel('c_{l,m}'); grid on;
end
