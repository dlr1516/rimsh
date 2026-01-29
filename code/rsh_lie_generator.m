function [Lx, Ly, Lz] = rsh_lie_generator(l)
  sz = 2*l+1;       % size of matrices: 2*l+1
  Lx = zeros(sz);
  Ly = zeros(sz);
  Lz = zeros(sz);

  % rows of secondary diagonal
  idx_u = sub2ind([sz, sz], 1:sz-1, sz-1:-1:1);
  idx_d = sub2ind([sz, sz], 2:sz, sz:-1:2);
  m = -l:l-1;
  signs = [-ones(1,l), ones(1,l)];
  u_lm = signs .* sqrt((l-m) .* (l+m+1)) ./ 2
  u_lm(l:l+1) = u_lm(l:l+1) .* sqrt(2);
  d_lm = u_lm;
  d_lm(l:l+1) = 0;
  Lx(idx_u) = u_lm;
  Lx(idx_d) = d_lm;
  disp(Lx)
  
  Lz(sub2ind([sz,sz], 1:sz, sz:-1:1)) = -l:l;
end