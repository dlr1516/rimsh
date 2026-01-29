function [Lx, Ly, Lz] = rsh_lie_generator(l)
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