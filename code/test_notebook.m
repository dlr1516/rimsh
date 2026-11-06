% Generators of representation in l-space of Lie algebra so(3)
lmax = 3;
[Lx, Ly, Lz] = rsh_lie_generator(lmax);

% Rotation A: Lie algebra representation
ta = pi/180 * (46);
pa = pi/180 * (84);
a = [cos(pa)*sin(ta), sin(pa)*sin(ta), cos(ta)]';
theta_a = (pi/180)*(53);
Ra = axang2rotm([a' theta_a]);
A = theta_a * (a(1) * Lx + a(2) * Ly + a(3) * Lz);
Ma = expm(A);
invMa = expm(-A);

% Rotation B: Lie algebra representation
tb = pi/180 * (-22);
pb = pi/180 * (224);
b = [cos(pb)*sin(tb), sin(pb)*sin(tb), cos(tb)]';
theta_b = (pi/180)*(-34);
Rb = axang2rotm([b' theta_b]);
B = theta_b *(b(1) * Lx + b(2) * Ly + b(3) * Lz);
Mb = expm(B);
invMb = expm(-B);

% Expression equivalence
sim_prod_AB_1 = Ma * B * invMa;
lb_AB = (A * B - B * A);
lb_AAB = (A * lb_AB - lb_AB * A);
sim_prod_AB_2 = B + sin(theta_a) / theta_a * lb_AB + (1-cos(theta_a)) / (theta_a^2) * lb_AAB;
sim_prod_AB_3 = expm(lb_AB);
disp('exp(B) * A * exp(-B):')
disp(sim_prod_AB_1)
disp('B + sin(a)/a * [A,B] + (1-cos(a))/a^2 * [A, [A,B]]:')
disp(sim_prod_AB_2)
disp('exp([A,B]):')
disp(sim_prod_AB_3)

% Finding C s.t. exp(C) A exp(-C) is anti-diagonal
z = [0, 0, 1]';
c = cross(a,z);
c = c / norm(c);
theta_c = acos(a(3));
C = theta_c * (c(1) * Lx + c(2) * Ly + c(3) * Lz);
Mc = expm(C);
A_diag = Mc * A * Mc' ./ theta_a;
disp('anti-diagonalized A, A_diag:')
disp(A_diag)

% Antidiagonalization using [F] = rsh_rot_to_z(l, u)
[Fa] = rsh_rot_to_z(lmax, a); 
A_diag2 = Fa' * A * Fa ./ theta_a;
disp('anti-diagonalized A: A_diag2')
disp(A_diag2)

% Direct computation of exp(A)
I = eye(size(A));
A2 = A*A;
A3 = A2*A;
Ma2 = I ...
     + (sin(theta_a)/theta_a) * A ...
     + ((1-cos(theta_a))/theta_a^2) * A2 ...
     + ((theta_a - sin(theta_a))/theta_a^3) * A3;
disp(Ma)
disp(Ma2)

% Verification addition theorem
theta1 = pi/180 * (30);
phi1 = pi/180 * (-123);
theta2 = pi/180 * (-42);
phi2 = pi/180 * (86);
u1 = [sin(theta1) * cos(phi1); sin(theta1) * sin(phi1); cos(theta1)];
u2 = [sin(theta2) * cos(phi2); sin(theta2) * sin(phi2); cos(theta2)];
Y1 = rsh_eval(lmax, u1);
Y2 = rsh_eval(lmax, u2);
lt = 2;
Y1_l = Y1(sh_lm_to_index(lt,-lt):sh_lm_to_index(lt,lt));
Y2_l = Y2(sh_lm_to_index(lt,-lt):sh_lm_to_index(lt,lt));
% P_l = legendre_poly(lt, u1' * u2);
sum_l = 4*pi / (2 * lt + 1) * (Y1_l' * Y2_l);
disp('u1^T * u2 | acos(u1^T * u2) [deg]')
disp([u1' * u2, 180/pi*acos(u1' * u2)])
% disp('P_l')
% disp(P_l)
disp('4*pi / (2 * lt + 1) * (Y1_l^T * Y2_l)')
disp(sum_l')




% Testing bounds of 
ars_filename_src = '/home/dario/workspace/rimsh/ars3d_2026_06_24_15_46_19/ars3d_cloud_dragon_vrip_04_noise000_occl000_rand000_20260605_1535_59.mat';
load(ars_filename_src, "ars_coeffs");
load(ars_filename_src, "R_gt");
ars_coeffs_src = ars_coeffs;
R_src = R_gt;

ars_filename_dst = '/home/dario/workspace/rimsh/ars3d_2026_06_24_15_46_19/ars3d_cloud_dragon_vrip_07_noise000_occl000_rand000_20260605_1535_59.mat';
load(ars_filename_dst, "ars_coeffs");
load(ars_filename_dst, "R_gt");
ars_coeffs_dst = ars_coeffs;
R_dst = R_gt;

