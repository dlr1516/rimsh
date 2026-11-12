function F3 = rsh_rot_to_z5(l, u)
  u = u(:); n = 2*l + 1;
  if l == 0, F3 = 1; return; end
  theta = norm(u);
  if theta < 1e-14, F3 = eye(n); return; end

  [Lx, Ly, Lz] = rsh_lie_generator(l);
  w = u / theta;                                  % unit axis
  Lu = w(1)*Lx + w(2)*Ly + w(3)*Lz;
  rho = hypot(w(1), w(2));

  if rho > 1e-14
    L1 = (-w(2)*Lx + w(1)*Ly) / rho;
    L2 = (-w(1)*w(3)*Lx - w(2)*w(3)*Ly + rho^2*Lz) / rho;
  else
    L1 = Lx;  L2 = sign(w(3)+eps) * Ly;
  end

  % v0: real null vector of Lu (use your rsh_eval-based a0 if norm(Lu*a0)~0)
  [~,~,VV] = svd(Lu);
  v = VV(:, end);  v = v / norm(v);
  v = fixphase(v);                                % real, largest component > 0

  % alternative computation
  % v_all = rsh_eval(l, u);
  % idx0 = sh_lm_to_index(l, -l:l);
  % a0 = v_all(idx0);
  % a0 = a0(:);
  % a0 = a0 / norm(a0);
  % v = a0

  F3 = zeros(n, n);
  F3(:, l+1) = v;                                  % v0 column (no sqrt(2))

  for m = 0:(l-1)
    Nm = sqrt((l-m)*(l+m+1));
    v = (L1 - 1i*L2) * v / Nm;                  % Lu v = i(m+1) v
    v = v / norm(v);                            % re-normalize
    vp = fixphase(v);                           % pin the phase like eig
    F3(:, l+1-(m+1)) = sqrt(2) * real(vp);       % a_{m+1}
    F3(:, l+1+(m+1)) = sqrt(2) * imag(vp);       % b_{m+1}
    % keep recurring from the un-rotated v (phase does not matter for the next step)
  end
end

function v = fixphase(v)
  [~, k] = max(abs(v)); 
  v = v * exp(-1i*angle(v(k)));                   % largest component real, positive
end