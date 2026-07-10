function ub = rsh_corr_ub_l(l, d, s, x_low, x_upp)
% Upper bound of f(a) = d' * M_l(R(a)) * s over a Rodrigues box.
%
% Inputs
%   l        : band index
%   d          : (2*l+1)x1 vector
%   s          : (2*l+1)x1 vector
%   x_low      : 3x1 box lower values
%   x_upp      : 3x1 box upper values

%
% Output
%   ub         : scalar upper bound on max_{a in box} d' * M_l(R(a)) * s

    x_mid = (x_upp + x_low) / 2;
    x_dim = (x_upp - x_low) / 2;
    rho_box = norm(x_dim);

    nd = norm(d);
    ns = norm(s);

    if nd == 0 || ns == 0
        ub = 0;
        return;
    end

    % Center value: u0 = M_l(R(a0)) s
    M_l = rsh_rot_lie(l, x_mid);
    t_mid = M_l * s;

    f_mid = d' * t_mid;
    beta_box = min(pi, l * rho_box);
    
    
    cos_mid = f0_mid / (nd * ns);   
    cos_mid = max(-1, min(1, cos_mid));
    theta_mid = acos(cos_mid);
    if theta0_mid <= beta_box
        ub = nd * ns;
    else
        rad = max(0, nd^2 * ns^2 - f0^2);
        ub = f_mid * cos(beta_box) + sqrt(rad) * sin(beta_box);
    end
    ub = min(ub, nd * ns);
end