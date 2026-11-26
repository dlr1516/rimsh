function c = von_mises_fisher_rsh_coeffs(lmax, mu, kappa)
%VON_MISES_FISHER_RSH_COEFFS  Real SH coefficients of a von Mises-Fisher density.
%
%   c = VON_MISES_FISHER_RSH_COEFFS(lmax, mu, kappa)
%
%   Closed-form coefficients (no numerical quadrature) of
%   f(u) = kappa/(4*pi*sinh(kappa)) * exp(kappa * mu'*u)
%   c_lm = int f(u) Y_lm(u) du = lambda_l(kappa) * Y_lm(mu)
%   lambda_l = kappa * i_l(kappa) / sinh(kappa)
%
%   c is a column vector of length (lmax+1)^2, packed in the same order
%   as rsh_eval (index of (l,m) = l^2 + l + m + 1).
%
%   See also RSH_EVAL, RSH_ROT_COEFFS.
  if kappa < 0, error('kappa must be >= 0'); end
  mu = mu(:) / norm(mu);

  % radial part (overflow-safe)
  lam = zeros(lmax+1, 1);
  if kappa == 0
    lam(1) = 1;
  else
    l   = (0:lmax).';
    Is  = besseli(l + 0.5, kappa, 1);              % I_{l+1/2}(kappa)*exp(-kappa)
    lam = sqrt(pi*kappa/2) * Is * 2 ./ (-expm1(-2*kappa));
  end

  % angles of mu: theta = polar angle (colatitude), phi = azimuth
  theta = acos(max(-1, min(1, mu(3))));
  phi   = atan2(mu(2), mu(1));

  Y = rsh_eval(lmax, theta, phi);                    % (lmax+1)^2 values
  Y = Y(:);

  % expand lambda_l to every (l,m) entry
  idx = zeros((lmax+1)^2, 1);
  for l = 0:lmax
    idx(l^2+1 : (l+1)^2) = l + 1; 
  end
  c = lam(idx) .* Y;
end
