function [Out] = rsh_rot_ivanic_v(m, n, l, R, Mprev)
%RSH_ROT_IVANIC_V  Function V of the Ivanic-Ruedenberg recursion.
%
%   Out = RSH_ROT_IVANIC_V(m, n, l, R, Mprev)
%
%   Helper of RSH_ROT_IVANIC (Ivanic & Ruedenberg, J. Phys. Chem. 1996/1998).
%   l is the band, R the 3x3 band-1 matrix, Mprev the matrix of band l-1.
% Compute two Kronecker deltas about m
if abs(m) == 1
    deltaM1 = 1;
else
    deltaM1 = 0;
end

if (m == 0)
    %  V = P(1, 1, n, l, r) + P(-1, -1, n, l, r);
    Out = rsh_rot_ivanic_p(1, 1, n, l, R, Mprev) + rsh_rot_ivanic_p(-1, -1, n, l, R, Mprev);
    return;
elseif (m > 0)
    %       P(1, m - 1, n, l, r) * sqrt(1 + KroneckerDelta(m, 1)) -
    %       P(-1, -m + 1, n, l, r) * (1 - KroneckerDelta(m, 1));
    Out = rsh_rot_ivanic_p(1, m - 1, n, l, R, Mprev) * sqrt(1 + deltaM1) - ...
          rsh_rot_ivanic_p(-1, -m + 1, n, l, R, Mprev) * (1 - deltaM1);
    return;
else
    %      P(1, m + 1, n, l, r) * (1 - KroneckerDelta(m, -1)) +
    %      P(-1, -m - 1, n, l, r) * sqrt(1 + KroneckerDelta(m, -1));
    Out = rsh_rot_ivanic_p(1, m + 1, n, l, R, Mprev) * (1 - deltaM1) + ... 
          rsh_rot_ivanic_p(-1, -m - 1, n, l, R, Mprev) * sqrt(1 + deltaM1);
end
end
