function [x_init, R_init] = rsh_corr_init_guess(lmax, rsh_src, rsh_dst)
% RSH_CORR_INIT_GUESS
% It estimates the value of the rotation that better aligns two functions
% on sphere expanded in Real Spherical Harmonics (RSH). 
% The estimates is meant as an initial guess for a local search optimization
% such as Gauss-Newton method. 
%
% MOTIVATION
% Let s and d be the 2*l+1 RSH coefficients of band l.
% The two coeffient vectors are related as d = M(R) * s where M(R) is the
% rotation matrix. 
% Using the representation of the Lie algebra so(3) generators Lx, Ly and Lz
% for band l, then 
%   M(R) = expm(x(1) * Lx + x(2) * Ly + x(3) * Lz).
% Expand and truncate the above exponential at first order as:
%   M(R) ~= I + (x(1) * Lx + x(2) * Ly + x(3) * Lz)
% Rewrite the equation d = M(R) * s as 
%   I * s + (x(1) * Lx + x(2) * Ly + x(3) * Lz) * s = d
%   [Lx*s, Ly*s, Lz*s] * x = d - s
% The initial guess is found by solving the linear system with (2*l+1) 
% equations in 3 unknown variables x.
%
% SUGGESTION
% The method operates on the coeffients of a single band l. 
% The user must select  
%
% INPUTS:
% - l: the band size of RSH to be used in optmization
% - rsh_src: RSH coefficients of source function
% - rsh_dst: RSH coefficients of destination function
%
% OUTPUTS:
% - x_init: Rodrigues paramaters of rotation, i.e. the axis vector multiplied 
%           by the rotation angle in radians;
% - R_init: the the initial guess in the form of rotation matrix. 
%
  A = zeros(0,3);
  b = [];
  for l=1:lmax
    [Lx, Ly, Lz] = rsh_lie_generator(l);
    s = rsh_src(sh_lm_to_index(l,-l:l));
    d = rsh_dst(sh_lm_to_index(l,-l:l));
    A = [A; 
         Lx * s, Ly * s, Lz * s];
    b = [b; 
         d - s];
  end
  x_init = A \ b;
  R_init = exp_SO3(x_init); 
return
