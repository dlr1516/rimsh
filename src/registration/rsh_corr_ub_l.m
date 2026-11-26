function [ub, lb] = rsh_corr_ub_l(l, s, d, x_low, x_upp)
%RSH_CORR_UB_L  Bounds of the band-l correlation over a box of Rodrigues vectors.
%
%   [ub, lb] = RSH_CORR_UB_L(l, s, d, x_low, x_upp)
%
%   Upper (and trivial lower) bound of f(a) = d' * M_l(R(a)) * s over the box
%   x_low <= a <= x_upp, e.g. for branch-and-bound search on SO(3).
%   NOTE: experimental code, not used in the experiments of the paper.
%
% Inputs
%   l          : band index
%   s          : (2*l+1)x1 vector of source RSH coefficient
%   d          : (2*l+1)x1 vector of destination RSH coefficient
%   x_low      : 3x1 box lower values
%   x_upp      : 3x1 box upper values
%
% Output
%   ub         : scalar upper bound on max_{a in box} d' * M_l(R(a)) * s
%   lb         : value at the box centre (a lower bound on the maximum)

    % column vectors
    s = s(:);  
    d = d(:); 
    x_low = x_low(:);
    x_upp = x_upp(:);

    % Box parameters
    x_mid = (x_upp + x_low) / 2;
    x_dim = (x_upp - x_low) / 2;
    rho_box = norm(x_dim);
    nd = norm(d);
    ns = norm(s);
    if abs(nd) < 1e-9 || abs(ns) < 1e-9
      ub = 0;
      lb = 0;
      return;
    end

    % Center value: u0 = M_l(R(a0)) s
    M_l = rsh_rot_lie(l, x_mid);
    t_mid = M_l * s;
    f_mid = d' * t_mid;

    % Maximum beta
    beta_box = min(pi, l * rho_box);    
    cos_mid = f_mid / (nd * ns);   
    cos_mid = max(-1, min(1, cos_mid));  % saturation on cos_mid (in principle not required)
    theta_mid = acos(cos_mid);
    if theta_mid <= beta_box
        ub = nd * ns;
    else
        rad = max(0, nd^2 * ns^2 - f_mid^2);
        ub = f_mid * cos(beta_box) + sqrt(rad) * sin(beta_box);
    end
    lb = f_mid;
    ub = min(ub, nd * ns);
end
