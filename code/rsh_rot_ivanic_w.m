function [Out] = rsh_rot_ivanic_w(m, n, l, R, Mprev)
if (m == 0)
    Out = 0;
elseif (m > 0)
    % P(1, m + 1, n, l, r) + P(-1, -m - 1, n, l, r);
    Out = rsh_rot_ivanic_p(1, m + 1, n, l, R, Mprev) + rsh_rot_ivanic_p(-1, -m - 1, n, l, R, Mprev);
else
    % P(1, m - 1, n, l, r) - P(-1, -m + 1, n, l, r);
    Out = rsh_rot_ivanic_p(1, m - 1, n, l, R, Mprev) - rsh_rot_ivanic_p(-1, -m + 1, n, l, R, Mprev);
end

end 
