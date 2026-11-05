%fclcunction test_rsh_lie_generator()
  lmax = 3;

  % Creates a rotation given its rotation axis and angle
  tw = pi/180 * (52);
  pw = pi/180 * (128);
  rot_axis_w = [cos(pw)*sin(tw), sin(pw)*sin(tw), cos(tw)]';
  %rot_angle_w = (pi/180) * (25);
  rot_angle_w = 1;
  Rw = axang2rotm([rot_axis_w' rot_angle_w]);

  % using Ivanic rotation formula
  [Ml_w1, M_w1] = rsh_rot_ivanic(Rw, lmax);
  disp('Ml_w according to Ivanic');
  disp(Ml_w1);

  % using Pinchon rotation formula
  % eul_zyz = rotm2eul(Rw, 'ZYZ');
  %[Ml_w3, M_w3] = rsh_rot_pinchon(eul_zyz, lmax);
  %disp('Ml_w according to Pinchon');
  %disp(Ml_w3);

  % Using Lie algebra
  [Lx, Ly, Lz] = rsh_lie_generator(lmax);
  Omega = (rot_axis_w(1) * Lx + rot_axis_w(2) * Ly + rot_axis_w(3) * Lz);
  %[V,J] = jordan(Omega)
  Ml_w2 = expm(rot_angle_w * Omega);
  disp('Ml_w as exponent of Lie algebra');
  disp(Ml_w2);

  % Using Lie algebra Eigenvectors
  Ml_w3 = rsh_rot_lie_eig(lmax, rot_angle_w * rot_axis_w);
  disp('Ml_w as exponent of Lie algebra Eigenvector decomposition');
  disp(Ml_w3);

  % Using Tchebychev
  %Ml_w3 = rsh_rot_tchebychev(lmax, rot_angle_w * rot_axis_w);
  %disp('Ml_w in Tchebychev expansion');
  %disp(Ml_w3);

  % Relation Omege + t^2 Omega = 0 (t short for rot_angle_w)
  % ra = rot_angle_w;
  % for k = 1:5
  %   disp(['[Omega | Omega^' num2str(k)])
  %   disp([Omega, Omega^k])
  % end
  