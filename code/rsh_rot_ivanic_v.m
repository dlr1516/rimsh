function [Out] = rsh_rot_ivanic_v(m, n, l, R, Mprev)
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
