function [x_init, R_init] = rsh_corr_init_guess(lmax, rsh_src, rsh_dst)
%RSH_CORR_INIT_GUESS  Linearized initial guess of the rotation aligning two RSH expansions.
%
%   [x_init, R_init] = RSH_CORR_INIT_GUESS(lmax, s, d)
%
% It estimates the value of the rotation that best aligns two functions
% on sphere expanded in Real Spherical Harmonics (RSH). 
% The estimate is meant as an initial guess for a local search optimization
% such as Gauss-Newton method. 
%
% MOTIVATION
% Let s and d be the 2*l+1 RSH coefficients of band l.
% The two coefficient vectors are related as d = M(R) * s where M(R) is the
% rotation matrix. 
% Using the representation of the Lie algebra so(3) generators Lx, Ly and Lz
% for band l, then 
%   M(R) = expm(x(1) * Lx + x(2) * Ly + x(3) * Lz).
% Expand and truncate the above exponential at first order as:
%   M(R) ~= I + (x(1) * Lx + x(2) * Ly + x(3) * Lz)
% Rewrite the equation d = M(R) * s as 
%   I * s + (x(1) * Lx + x(2) * Ly + x(3) * Lz) * s = d
%   [Lx*s, Ly*s, Lz*s] * x = d - s
% The initial guess is the least-squares solution of the stacked linear
% systems of bands l = 1..lmax (3 unknowns). The first-order expansion is only
% accurate for small rotations.
%
% INPUTS:
% - lmax: maximum band of RSH used in the estimation
% - rsh_src: RSH coefficients of source function
% - rsh_dst: RSH coefficients of destination function
%
% OUTPUTS:
% - x_init: Rodrigues parameters of rotation, i.e. the axis vector multiplied 
%           by the rotation angle in radians;
% - R_init: the initial guess in the form of rotation matrix. 
%
  A = zeros(0,3);
  b = [];
  for l=1:lmax
    [Lx, Ly, Lz] = rsh_lie_generator(l);
    s = reshape(rsh_src(sh_lm_to_index(l,-l:l)), [], 1);
    d = reshape(rsh_dst(sh_lm_to_index(l,-l:l)), [], 1);
    A = [A; 
         Lx * s, Ly * s, Lz * s];
    b = [b; 
         d - s];
  end
  x_init = A \ b;
  R_init = exp_SO3(x_init);
end
