function F = rsh_rot_to_z(l, u, G)
%RSH_ROT_TO_Z  Real orthogonal factor of the anti-diagonal factorization of RSH rotation.
%
%   F = rsh_rot_to_z(l, u)
%   F = rsh_rot_to_z(l, u, G)
%
% Inputs
%   l : band (degree)
%   u : 3-vector, Rodrigues vector or rotation axis (only its direction is used)
%   G : optional struct with fields Lx, Ly, Lz (generators of band l, as returned by
%       rsh_lie_generator(l)). Pass it when calling repeatedly: building the generators
%       is otherwise repeated at every call.
%
% Output
%   F : (2l+1)x(2l+1) orthogonal matrix such that, for u = theta*u_hat,
%           L_{l,u_hat} = F * L_z * F'      and      M_l(u) = F * M_lz(theta) * F'
%       Column order (indices -l..l):
%           [sqrt(2) a_l, ..., sqrt(2) a_1, a_0, sqrt(2) b_1, ..., sqrt(2) b_l]
%       where v_m = a_m + 1i*b_m is the eigenvector of L_{l,u_hat} with eigenvalue 1i*m.
%
% Algorithm
%   v_0 = Y_l(u_hat) (normalized); v_{m+1} = (L2 + 1i*L1) v_m / gamma_{l,m}, with
%   L1 = e1.L, L2 = e2.L, (e1, e2, u_hat) right-handed and gamma_{l,m} = sqrt((l-m)(l+m+1)).
%   The plain recurrence amplifies round-off by roughly prod gamma_{l,k}/gamma_{l,m}
%   (||F'F - I|| ~ 1e-11 at l = 32, 1e-8 at l = 64). Each new pair (a_m, b_m) is
%   therefore re-orthogonalized against all previous ones (classical Gram-Schmidt, two
%   passes), which keeps F orthogonal to round-off level at O(l^3) cost per band.
%
% Requires: rsh_eval, sh_lm_to_index, rsh_lie_generator (if G is not given).
% See also RSH_ROT_LIE_EIG.
% Note: the accuracy of F is limited by the accuracy of Y_l(u_hat) from rsh_eval. If
% rsh_eval computes sin(theta) as sqrt(1 - cos(theta)^2), axes close to +/-z lose
% accuracy (error ~ eps/rho); compute sin(theta) = hypot(u_x, u_y)/norm(u) instead.

    u = u(:);
    if numel(u) ~= 3, error('u must be a 3-vector.'); end
    n = 2*l + 1;
    if l == 0, F = 1; return; end
    theta = norm(u);
    if theta < 1e-14, F = eye(n); return; end      % no axis: any F works, M = I
    if nargin < 3 || isempty(G)
        [Lx, Ly, Lz] = rsh_lie_generator(l);
    else
        Lx = G.Lx; Ly = G.Ly; Lz = G.Lz;
    end
    ux = u(1); uy = u(2); uz = u(3);
    rho = hypot(ux, uy);

    % L1 = e1.L, L2 = e2.L with e1 = (-uy, ux, 0)/rho, e2 = (-ux uz, -uy uz, rho^2)/(theta rho)
    % (both unit vectors for the NON-normalized u; e1 x e2 = u/theta).
    % The generators have at most 2 nonzeros per row: sparse products cost O(l).
    if rho > 1e-14
        L1 = sparse((-uy * Lx + ux * Ly) / rho);
        L2 = sparse((-ux * uz * Lx - uy * uz * Ly + rho^2 * Lz) / (theta * rho));
    else                                            % axis = +/- z
        L1 = sparse(Lx);
        if uz >= 0, L2 = sparse(Ly); else, L2 = -sparse(Ly); end
    end

    % v_0 = Y_l(u_hat)
    v_all = rsh_eval(l, u);
    a0 = v_all(sh_lm_to_index(l, -l:l));
    a0 = a0(:) / norm(a0);

    % Q: orthonormal basis of span{a_j, b_j} computed so far (columns 1..k)
    Q = zeros(n, n); Q(:,1) = a0; k = 1;
    A = zeros(n, l); B = zeros(n, l);               % A(:,m) = a_m, B(:,m) = b_m
    a = a0; b = zeros(n, 1);
    for m = 0:(l-1)
        Nm = sqrt((l - m) * (l + m + 1));
        an = (L2 * a - L1 * b) / Nm;
        bn = (L2 * b + L1 * a) / Nm;
        Qk = Q(:, 1:k);
        for pass = 1:2                              % CGS2 re-orthogonalization
            an = an - Qk * (Qk' * an);
            bn = bn - Qk * (Qk' * bn);
        end
        s = sqrt(an' * an + bn' * bn);              % ||v_{m+1}|| = 1
        if s < 1e-10, error('Degenerate recurrence at m = %d.', m+1); end
        a = an / s; b = bn / s;
        A(:, m+1) = a; B(:, m+1) = b;
        Q(:, k+1) = a / norm(a); Q(:, k+2) = b / norm(b); k = k + 2;
    end

    F = zeros(n, n);
    F(:, l+1) = a0;
    F(:, l:-1:1)   = sqrt(2) * A;                   % column l+1-m <- sqrt(2) a_m
    F(:, l+2:end)  = sqrt(2) * B;                   % column l+1+m <- sqrt(2) b_m
end
