function [f, g, H] = rsh_corr_grad_hess(lmax, rsh_src, rsh_dst, x)
%RSH_CORR_GRAD_HESS  Correlation of two RSH expansions, gradient and Hessian on SO(3).
%
%   [f, g, H] = RSH_CORR_GRAD_HESS(lmax, s, d, x)
%
%   Evaluates the correlation between a source and a destination function on
%   the sphere, given by their RSH coefficients s and d,
%       f(R) = sum_{l=1..lmax} d_l' * M_l(R) * s_l,
%   together with its closed-form gradient and Hessian with respect to a right
%   perturbation R * exp([w]_x) at w = 0:
%       g_k  = sum_l dt_l' * L_k * s_l,
%       H_jk = sum_l dt_l' * (L_j L_k + L_k L_j) / 2 * s_l,   dt_l = M_l(R)' * d_l.
%   Band l = 0 is rotation invariant and is not included in f.
%
%   Inputs
%     lmax : maximum band
%     s, d : (lmax+1)^2 vectors of source / destination RSH coefficients
%     x    : rotation, either a 3x3 matrix or a 3-vector Rodrigues vector
%
%   Outputs
%     f : correlation value
%     g : 3x1 gradient
%     H : 3x3 (symmetric) Hessian
%
%   See also RSH_OPT_NEWTON_SO3, RSH_LIE_GENERATOR.

  if isequal(size(x), [3 3])
    x = log_SO3(x);
  elseif numel(x) == 3
    x = x(:);
  else
    error('rsh_corr_grad_hess:input', ...
          'x must be a 3x3 rotation matrix or a 3-vector of Lie algebra parameters.');
  end

  f = 0;
  g = zeros(3, 1);
  H = zeros(3, 3);

  for l = 1:lmax
    L = zeros(2*l+1, 2*l+1, 3);
    [L(:,:,1), L(:,:,2), L(:,:,3)] = rsh_lie_generator(l);
    % coefficients of band l
    idx = sh_lm_to_index(l, -l:l);
    s = reshape(rsh_src(idx), [], 1);
    d = reshape(rsh_dst(idx), [], 1);
    M = expm(x(1) * L(:,:,1) + x(2) * L(:,:,2) + x(3) * L(:,:,3));
    dt = M' * d;
    f = f + d' * M * s;
    % gradient
    for j = 1:3
      g(j) = g(j) + dt' * L(:,:,j) * s;
    end
    % Hessian (symmetrized second derivative, diagonal included)
    for j = 1:3
      for k = j:3
        Inc = dt' * (L(:,:,j) * L(:,:,k) + L(:,:,k) * L(:,:,j)) * s / 2;
        H(j,k) = H(j,k) + Inc;
        if j ~= k
          H(k,j) = H(k,j) + Inc;
        end
      end
    end
  end
end
