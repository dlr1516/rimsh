function test_rsh_corr_grad_hess()
  % reading and setting data
  lmax = 6;
  disp('reading src coeffs:')
  % file_src = fopen('../data/ars_coeffs_src.bin','r');
  % if file_src == -1
  %   error('Failed to open src file');
  % end
  % coeffs_src = fread(file_src);
  % fclose(file_src);
  file_src = matfile('../data/ars_coeffs_src.mat');
  coeffs_src = file_src.ars_src;

  disp('reading dst coeffs:')
  % file_dst = fopen('../data/ars_coeffs_dst.bin','r');
  % if file_dst == -1
  %   error('Failed to open dst file');
  % end
  % coeffs_dst = fread(file_dst);
  % fclose(file_dst);
  file_dst = matfile('../data/ars_coeffs_dst.mat');
  coeffs_dst = file_dst.ars_dst;

  % Groundtruth:
  tw = pi/180 * (52);
  pw = pi/180 * (128);
  rot_axis_gt = [cos(pw)*sin(tw), sin(pw)*sin(tw), cos(tw)]';
  rot_angle_gt = (pi/180) * (25);
  Rgt = axang2rotm([rot_axis_gt' rot_angle_gt]);
  coeffs_rot = zeros(length(coeffs_dst),1);
  % Test groundtruth
  [~, Mgt] = rsh_rot_ivanic(Rgt, lmax);
  for l = 1:lmax
    idx = sh_lm_to_index(l,-l:l);
    s = coeffs_src(idx);
    d = coeffs_dst(idx);
    d_pred = Mgt{l+1} * s;
    coeffs_rot(idx) = d_pred;
    disp(['l ' num2str(l) ': coeffs src | rotated src | dst'])
    disp([s, d_pred, d]);

    %[Lx, Ly, Lz] = rsh_lie_generator(l);
    %Omega = (rot_axis_gt(1) * Lx + rot_axis_gt(2) * Ly + rot_axis_gt(3) * Lz);
    %Mlie = expm(rot_angle_gt * Omega);
    %disp('Mgt | Mlie')
    %disp([Mgt{l+1}; Mlie])
  end

  % Computation of gradient
  [g,H] = rsh_corr_grad_hess(coeffs_src, coeffs_dst, lmax);
  disp('gradient')
  disp(g)
  disp('Hessian')
  disp(H)
  disp('eigenvalues of Hessian')
  disp(eig(H))

  %
  disp('gt rotation lie params:')
  disp(rot_angle_gt * rot_axis_gt)
  x = H \ (-g)

end