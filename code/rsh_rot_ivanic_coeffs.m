function [u, v, w] = rsh_rot_ivanic_coeffs(m, n, l)
  % d = KroneckerDelta(m, 0);
  % Compute two Kronecker deltas about m
  if m == 0
    deltaM0 = 1;
  else
    deltaM0 = 0;
  end


  % denom = (abs(n) == l ? 2.0 * l * (2.0 * l - 1) : (l + n) * (l - n));
  if abs(n) == l 
      denom = 2 * l * (2.0 * l - 1);
  else
      denom =  (l + n) * (l - n);
  end

  % u = sqrt((l + m) * (l - m) / denom);
  u = sqrt((l + m) * (l - m) / denom);

  % v = 0.5 * sqrt((1 + d) * (l + abs(m) - 1.0) * (l + abs(m)) / denom) * (1 - 2 * d);
  v = 0.5 * sqrt((1 + deltaM0) * (l + abs(m) - 1.0) * (l + abs(m)) / denom) * (1 - 2 * deltaM0);

  % w = -0.5 * sqrt((l - abs(m) - 1) * (l - abs(m)) / denom) * (1 - d);
  w =-0.5 * sqrt((l - abs(m) - 1) * (l - abs(m)) / denom) * (1 - deltaM0);

end
