function [rsh_coeffs_rot] = rsh_rot_coeffs(R, rsh_coeffs)
  % Checks that the length of rsh_coeffs is equal to (lmax+1)^2
  len = length(rsh_coeffs);
  lmax = floor(sqrt(len)-1);
  assert(len == (lmax+1)^2)

  % 
  if size(R,1) == 3 & size(R,2) == 1
    x = R;
  else if size(R, 1) == 3 & size(R, 2) == 3
    x = log_SO3(R);
  else 
    error('input rotation R must be either a rotation matrix or a Rodrigues vector')
  end

  [~, M] = rsh_rot_lie(lmax, x);
  rsh_coeffs_rot = zeros(len,1);
  for l = 1:lmax
    idx = sh_lm_to_index(l,-l:l);
    c_l = rsh_coeffs(idx);
    c_l_rot = M{l+1} * c_l;
    rsh_coeffs_rot(idx) = c_l_rot;
  end
end