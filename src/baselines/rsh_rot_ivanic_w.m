function [Out] = rsh_rot_ivanic_w(m, n, l, R, Mprev)
%RSH_ROT_IVANIC_W  Function W of the Ivanic-Ruedenberg recursion.
%
%   Out = RSH_ROT_IVANIC_W(m, n, l, R, Mprev)
%
%   Helper of RSH_ROT_IVANIC (Ivanic & Ruedenberg, J. Phys. Chem. 1996/1998).
%   l is the band, R the 3x3 band-1 matrix, Mprev the matrix of band l-1.
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
