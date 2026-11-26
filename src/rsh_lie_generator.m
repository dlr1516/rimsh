function [Lx, Ly, Lz] = rsh_lie_generator(l)
%RSH_LIE_GENERATOR  Generators of so(3) in the real spherical harmonics basis of band l.
%
%   [Lx, Ly, Lz] = RSH_LIE_GENERATOR(l)
%
%   Returns the (2l+1)x(2l+1) real skew-symmetric matrices that represent the
%   generators of the Lie algebra so(3) acting on the real spherical harmonics
%   (RSH) coefficients of band l, ordered as m = -l, ..., l (same order as
%   RSH_EVAL and SH_LM_TO_INDEX).
%
%   With x = theta * u (Rodrigues / rotation vector) the RSH rotation matrix of
%   band l is
%       M_l(x) = expm(x(1)*Lx + x(2)*Ly + x(3)*Lz),
%   and it satisfies Y_l(R*u) = M_l(x) * Y_l(u) with R = exp([x]_x).
%   The generators fulfil [Lx, Ly] = Lz (and cyclic permutations) and
%   Lx^2 + Ly^2 + Lz^2 = -l(l+1) I. Each matrix has at most two non-zero
%   entries per row (Lx, Ly near the anti-diagonal, Lz on the anti-diagonal).
%
%   Example
%       [Lx, Ly, Lz] = rsh_lie_generator(2);
%       M2 = expm(0.3*Lx - 0.1*Ly + 0.5*Lz);   % band-2 rotation matrix
%
%   See also RSH_ROT_LIE_EIG, RSH_ROT_LIE, RSH_ROT_TO_Z.
  sz = 2*l+1;       % size of matrices: 2*l+1
  Lx = zeros(sz);
  Ly = zeros(sz);
  Lz = zeros(sz);
  m = -l:l-1;
  signs = [-ones(1,l), ones(1,l)];

  % Lx: 
  idx_u = sub2ind([sz, sz], 1:sz-1, sz-1:-1:1);
  idx_d = sub2ind([sz, sz], 2:sz, sz:-1:2);
  u_lm = signs .* sqrt((l-m) .* (l+m+1)) ./ 2;
  u_lm(l:l+1) = u_lm(l:l+1) .* sqrt(2);
  d_lm = u_lm;
  d_lm(l:l+1) = 0;
  Lx(idx_u) = u_lm;
  Lx(idx_d) = d_lm;

  % Ly: 
  idx_u = sub2ind([sz, sz], 1:sz-1, 2:sz);
  idx_d = sub2ind([sz, sz], 2:sz, 1:sz-1);
  d_lm = signs .* sqrt((l-m) .* (l+m+1)) ./ 2;
  d_lm(l) = 0;
  d_lm(l+1) = d_lm(l+1) * sqrt(2);
  u_lm = d_lm(end:-1:1);
  u_lm(l+1) = -u_lm(l);
  u_lm(l) = 0;
  Ly(idx_d) = d_lm;
  Ly(idx_u) = u_lm;
  
  % Lz: 
  Lz(sub2ind([sz,sz], 1:sz, sz:-1:1)) = l:-1:-l;
end
