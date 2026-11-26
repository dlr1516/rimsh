function rsh_coeffs_rot = rsh_rot_coeffs(R, rsh_coeffs)
%RSH_ROT_COEFFS  Rotate a real spherical harmonics expansion.
%
%   c_rot = RSH_ROT_COEFFS(R, c)
%
%   Given the RSH coefficients c of a function f on the sphere (packed as in
%   SH_LM_TO_INDEX, length (lmax+1)^2), returns the coefficients of the rotated
%   function f_rot(u) = f(R' * u), i.e. c_rot_l = M_l(R) * c_l for every band.
%
%   Inputs
%     R : 3x3 rotation matrix, or 3-vector Rodrigues / rotation vector
%     c : (lmax+1)^2 vector of RSH coefficients
%
%   Output
%     c_rot : rotated coefficients, same size as c
%
%   Example
%       c  = von_mises_fisher_rsh_coeffs(8, [0;0;1], 20);
%       cr = rsh_rot_coeffs([pi/2; 0; 0], c);   % vMF now centred on [0;-1;0]
%
%   See also RSH_ROT_LIE_EIG, RSH_EVAL.

  len = numel(rsh_coeffs);
  lmax = round(sqrt(len)) - 1;
  if len ~= (lmax+1)^2
    error('rsh_rot_coeffs:size', ...
          'length of rsh_coeffs (%d) must be (lmax+1)^2.', len);
  end

  if isequal(size(R), [3 3])
    x = log_SO3(R);
  elseif numel(R) == 3
    x = R(:);
  else
    error('rsh_rot_coeffs:input', ...
          'R must be a 3x3 rotation matrix or a 3-vector Rodrigues vector.');
  end

  [~, M] = rsh_rot_lie_eig(lmax, x);
  rsh_coeffs_rot = zeros(size(rsh_coeffs));
  rsh_coeffs_rot(1) = rsh_coeffs(1);              % band 0 is invariant
  for l = 1:lmax
    idx = sh_lm_to_index(l, -l:l);
    rsh_coeffs_rot(idx) = M{l+1} * reshape(rsh_coeffs(idx), [], 1);
  end
end
