function [F] = rsh_rot_to_z(l, u)
%RSH_ROT_TO_Z  Real orthogonal change of basis for RSH rotation blocks.
%
%   F = rsh_rot_to_z(l, u)
%
% Inputs
%   l : band/order degree
%   u : 3x1 or 1x3 Rodrigues rotation vector
%
% Output
%   F : (2*l+1)x(2*l+1) real matrix such that
%           M(u) = F * M_z * F.'
%       where M_z is the real block rotation matrix for rotations about z:
%           blocks R(m*theta), m = l,l-1,...,1 and scalar 1 for m=0.
%
% Notes
%   - Requires:
%       rsh_eval(...)
%       rsh_lm_to_index(...)
%       and generators Lx,Ly,Lz for the band l.
%   - Here it is assumed that a helper function
%         [Lx, Ly, Lz] = rsh_generators(l)
%     exists and returns the (2*l+1)x(2*l+1) real generators in the RSH basis.
%   - The column ordering of F is:
%         [sqrt(2)*a_l, sqrt(2)*b_l, sqrt(2)*a_{l-1}, sqrt(2)*b_{l-1}, ...
%          sqrt(2)*a_1, sqrt(2)*b_1, a_0]
%
%   with v_m = a_m + i b_m and U v_m = i m theta v_m.
    u = u(:);
    n = 2*l + 1;

    if l == 0
        Fa = 1; return;
    end

    theta = norm(u);
    if theta < 1e-14
        Fa = eye(n); return;
    end

    axis = u / theta;
    [Lx, Ly, Lz] = rsh_lie_generator(l);
    G = axis(1)*Lx + axis(2)*Ly + axis(3)*Lz;   % U / theta

    % v0
    v_all = rsh_eval(l, axis);
    idx   = rsh_lm_to_index(l, -l:l);
    a0    = v_all(idx);
    a0    = a0(:) / norm(a0);

    % complex eigenvectors of G
    [V, D] = eig(G);
    evals  = diag(D);

    a = cell(l,1);
    b = cell(l,1);

    for m = 1:l
        [~, j] = min(abs(evals - 1i*m));
        v = V(:,j);
        v = v / norm(v);

        a_m = real(v);
        b_m = imag(v);

        % orthonormalize the pair
        a_m = a_m / norm(a_m);
        b_m = b_m - (a_m' * b_m) * a_m;
        b_m = b_m / norm(b_m);

        % ---- sign rule: enforce G*a_m = -m*b_m ----
        if b_m' * G * a_m > 0
            b_m = -b_m;
        end

        a{m} = a_m;
        b{m} = b_m;
    end

    % ---- assembly: [b_l ... b_1, a_0, a_1 ... a_l] ----
    Fa  = zeros(n, n);
    col = 1;

    for m = l:-1:1
        Fa(:,col) = b{m};  col = col + 1;
    end
    Fa(:,col) = a0;        col = col + 1;
    for m = 1:l
        Fa(:,col) = a{m};  col = col + 1;
    end

    Fa = sqrt(2) * Fa;      % a_m, b_m already normalized individually
    Fa(:, l+1) = a0;        % a_0 must NOT be scaled

    % ---- final safety sign check on the anti-diagonal ----
    L = Fa' * G * Fa;
    for m = 1:l
        r = l + 1 - m;
        c = l + 1 + m;
        if L(r, c) < 0
            Fa(:, r) = -Fa(:, r);
        end
    end
end

