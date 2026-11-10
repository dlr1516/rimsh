function [F, V, D] = eig_to_basis(L)
  [n, ~] = size(L);
  l = floor((n-1)/2);

  [V, D] = eig(L);     % eigenvalues on diag(D), eigenvectors in V
  ev = diag(D);

  [~, idx] = sort(imag(ev), 'ascend');   % sort key: imaginary part
  D = D(idx, idx);
  V = V(:, idx);
  D = diag(diag(D));

  idx_up = (l+2):(2*l+1);
  idx_down = (2*l+1):-1:(l+2);
  c = real(V(:, (l+1)));
  A = sqrt(2) * real(V(:, idx_down));
  B = sqrt(2) * imag(V(:, idx_up));

  F = [A, c, B];             % (2n) x n : real parts over imag parts
end