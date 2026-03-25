% Generators of representation in l-space of Lie algebra so(3)
lmax = 3;
[Lx, Ly, Lz] = rsh_lie_generator(lmax);

% Rotation A: Lie algebra representation
ta = pi/180 * (52);
pa = pi/180 * (128);
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
B =theta_b *(b(1) * Lx + b(2) * Ly + b(3) * Lz);
Mb = expm(B);
invMb = expm(-B);

% Expression equivalence
sim_prod_AB_1 = Ma * B * invMa;
lb_AB = (A * B - B * A);
lb_AAB = (A * lb_AB - lb_AB * A);
sim_prod_AB_2 = B + sin(theta_a) / theta_a * lb_AB + (1-cos(theta_a)) / (theta_a^2) * lb_AAB;
disp('exp(B) * A * exp(-B):')
disp(sim_prod_AB_1)
disp('B + sin(a)/a * [A,B] + (1-cos(a))/a^2 * [A, [A,B]]:')
disp(sim_prod_AB_2)
