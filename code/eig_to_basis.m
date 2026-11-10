function [F, D, V] = eig_to_basis(L)
  [n, ~] = size(L);
  l = floor((n-1)/2)

  [V, D] = eig(L);     % eigenvalues on diag(D), eigenvectors in V
  ev = diag(D);

  [~, idx] = sort(imag(ev), 'ascend');   % sort key: imaginary part
  D = D(idx, idx);
  V = V(:, idx)
  D = diag(diag(D));

  c = real(V(:, (l+1)))
  A = sqrt(2) * real(V(:, (l+2):n))
  B = sqrt(2) * imag(V(:, 1:l))

  F = [B, c, A];             % (2n) x n : real parts over imag parts
end