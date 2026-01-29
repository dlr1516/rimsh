function test_rsh_rot_lie_algebra()
  lmax = 5;
  rot_angle = 0.1;
  rot_axis_x = [1 0 0]';
  Rx = axang2rotm([rot_axis_x' rot_angle]);
  [~, Ml_x] = rsh_rot_ivanic(Rx, lmax);
  Ysh_x = rsh_eval(lmax, rot_axis_x);
  %
  rot_axis_y = [0 1 0]';
  Ry = axang2rotm([rot_axis_y' rot_angle]);
  [~, Ml_y] = rsh_rot_ivanic(Ry, lmax);
  Ysh_y = rsh_eval(lmax, rot_axis_y);
  %
  rot_axis_z = [0 0 1]';
  Rz = axang2rotm([rot_axis_z' rot_angle]);
  [~, Ml_z] = rsh_rot_ivanic(Rz, lmax);
  Ysh_z = rsh_eval(lmax, rot_axis_z);
  for l = 0:lmax
      Yl_x = Ysh_x(sh_lm_to_index(l,-l):sh_lm_to_index(l,l));
      norm_Yl_x = norm(Yl_x);
      Yl_x = round(100000 * Yl_x ./ norm_Yl_x) / 100000;
      Yl_y = Ysh_y(sh_lm_to_index(l,-l):sh_lm_to_index(l,l));
      norm_Yl_y = norm(Yl_y);
      Yl_y = round(100000 * Yl_y ./ norm(Yl_y)) / 100000;
      Yl_z = Ysh_z(sh_lm_to_index(l,-l):sh_lm_to_index(l,l));
      norm_Yl_z = norm(Yl_z);
      Yl_z = round(100000 * Yl_z ./ norm(Yl_z)) / 100000;
      Ll_x = logm(Ml_x{l+1}) ./ rot_angle;
      Ll_y = logm(Ml_y{l+1}) ./ rot_angle;
      Ll_z = logm(Ml_z{l+1}) ./ rot_angle;
      [diag_x, diag_y, diag_z, diag_x2, diag_y2] = rsh_diagonals(Ll_x);
      disp(['********************** l=' num2str(l)])
      disp(['  Lx ' num2str(l) '  angle ' num2str(rot_angle) '  rot_axis [' num2str(rot_axis_x') ']'])
      disp(['      norm   RSh [' num2str(Yl_x') ']'])
      disp(['      unnorm Rsh ['  num2str(norm_Yl_x * Yl_x') ']'])
      disp(Ll_x)
      disp('diagonal x:')
      disp(diag_x)
      m = 0:l;
      u_lm = sqrt((l-m) .* (l+m+1)) ./ 2;
      u_lm(1) = u_lm(1) * sqrt(2);
      disp('u_lm (up, raising)')
      disp(u_lm);
      %[Vl_x,Dl_x] = eig_sort_imag(Ll_x);
      %disp('eigenvalues:') 
      %disp(diag(Dl_x)')
      disp('---')
      disp(['  Ly ' num2str(l) '  angle ' num2str(rot_angle) '  rot_axis [' num2str(rot_axis_y') ']'])
      disp(['      norm   RSh [' num2str(Yl_y') ']'])
      disp(['      unnorm Rsh ['  num2str(norm_Yl_y * Yl_y') ']'])
      disp('---')
      disp(Ll_y) 
      disp('diagonal y:')
      disp(diag_y)
      disp(u_lm);
      disp(u_lm./2)
      [Vl_y,Dl_y] = eig(Ll_y);
      %disp('eigenvalues:') 
      %disp(diag(Dl_y)')
      disp('---')
      disp(['  Lz ' num2str(l) '  angle ' num2str(rot_angle) '  rot_axis [' num2str(rot_axis_z') ']'])
      disp(['      clcnorm   RSh [' num2str(Yl_z') ']'])
      disp(['      unnorm Rsh ['  num2str(norm_Yl_z * Yl_z') ']'])
      disp(Ll_z)
      [Vl_z,Dl_z] = eig_sort_imag(Ll_z);
      %disp('eigenvalues:') 
      %disp(diag(Dl_z)')
  end

function [diag_x, diag_y, diag_z, diag_x2, diag_y2] = rsh_diagonals(M)
   s = size(M);
   if s(1) ~= s(2) || mod(s(1), 2) ~= 1
       error(['invalid size ' num2str(s) ': it must be square and odd-sized'])
   end
   n = s(1);
   l = (n-1)/2;
   idx = sub2ind(s, n:-1:1, 1:n);
   diag_z = M(idx);
   idx = sub2ind(s, 1:n-1, n-1:-1:1);
   diag_x = M(idx);
   idx = sub2ind(s, 2:n, n:-1:2);
   diag_x2 = M(idx);
   idx = sub2ind(s, 2:n, 1:n-1);
   diag_y = M(idx);
   idx = sub2ind(s, n-1:-1:1, n:-1:2);
   diag_y2 = M(idx);





