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
  Rgt = axang2rotm([rot_axis_gt', rot_angle_gt]);
  x_gt = log_SO3(Rgt);
  theta_gt = norm(x_gt);
  x_gt_antip = (pi - theta_gt) * x_gt / theta_gt;
  R_gt_antip = exp_SO3(x_gt_antip);

  coeffs_rot = zeros(length(coeffs_dst),1);
  % Test groundtruth
  [~, Mgt] = rsh_rot_ivanic(Rgt, lmax);
  [~, Mgt_antip] = rsh_rot_ivanic(R_gt_antip, lmax);
  for l = 1:lmax
    idx = sh_lm_to_index(l,-l:l);
    s = coeffs_src(idx);
    d = coeffs_dst(idx);
    d_pred = Mgt{l+1} * s;
    d_antip = Mgt_antip{l+1} * s;
    coeffs_rot(idx) = d_pred;
    disp(['l ' num2str(l) ': coeffs src | rotated src | dst | d_antip'])
    disp([s, d_pred, d, d_antip]);
    disp(['norm s ' num2str(norm(s)) ', d ' num2str(norm(d))])
  end

  disp('------------------------------------------------------------------')
  disp('Estimation of initial guess')
  [x_init, R_init] = rsh_corr_init_guess_2nd(lmax, coeffs_src, coeffs_dst);

  disp('Rgt')
  disp(Rgt)
  disp(['x_gt: angle[deg] ' num2str(180/pi * norm(x_gt)) ' axis'])
  disp(x_gt / norm(x_gt)) 


  disp('------------------------------------------------------------------')
  disp('Computing MAX of Correlation: R = max sum_{l} d_l^T * M_l(R) * s_l')
  
  % The initial guess x_init can be expressed in two forms (bot supported):
  % a) x_init as a 3D vector with the Lie parameters of rotation (aka 
  %    Rodrigues parameters, i.e. rotation axis multiplied by the angle);
  % b) x_init as a rotation matrix.
  % The user can set the most convenient one. 
  x_init = [0, 0, 0]';
  % Initial guess close to the known solution Rgt
  %ax_diff = [3 1 2];
  %ax_diff = ax_diff / norm(ax_diff);
  %x_init = Rgt * axang2rotm([ax_diff, (pi/180)*(5)])
  data.lmax = lmax;
  [Ropt, data] = rsh_opt_newton_so3(lmax, coeffs_src, coeffs_dst, x_init, data);
  disp('Ropt | Rgt | R_init')
  disp([Ropt, Rgt, R_init])
  [f_gt, g_gt, H_gt] = rsh_corr_grad_hess(lmax, coeffs_src, coeffs_dst, Rgt);
  disp(['function in Rgt: ' num2str(f_gt)])
  disp('gradient in Rgt, g_gt')
  disp(g_gt)

  rsh_plot(coeffs_src, 180, 'RSH src');
  figure
  rsh_plot(coeffs_dst, 180, 'RSH dst');
end