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
    if numel(u) ~= 3
        error('u must be a 3-vector.');
    end

    % Band dimension
    n = 2*l + 1;

    % Trivial case l = 0
    if l == 0
        F = 1;
        return;
    end

    % Rotation magnitude
    theta = norm(u);

    % No rotation: F is not unique. Return identity.
    if theta < 1e-14
        F = eye(n);
        return;
    end

    ux = u(1);
    uy = u(2);
    uz = u(3);
    rho = hypot(ux, uy);

    % RSH generators for this band
    [Lx, Ly, Lz] = rsh_lie_generator(l);

    Lu = (ux * Lx + uy * Ly + uz * Lz) / theta;

    % Build A and B defining L_+^(n) = A + i B
    if rho > 1e-14
        A = (-uy * Lx + ux * Ly) / rho;
        B = (-ux * uz * Lx - uy * uz * Ly + rho^2 * Lz) / (theta * rho);
    else
        % Axis = +/- z
        A = Lx;
        if uz >= 0
            B = Ly;
        else
            B = -Ly;
        end
    end

    % a0 = v0 from RSH evaluation along axis u
    v_all = rsh_eval(l, u);
    idx0 = sh_lm_to_index(l, -l:l);
    a0 = v_all(idx0);
    a0 = a0(:);

    % Normalize a0
    a0 = a0 / norm(a0);

    % Storage for a_m, b_m
    a = cell(l+1, 1);
    b = cell(l+1, 1);

    a{1} = a0;              % corresponds to m = 0
    b{1} = zeros(n, 1);     % b_0 = 0

    % Recurrence for m = 0,1,...,l-1
    for m = 0:(l-1)
        Nm = sqrt((l - m) * (l + m + 1));

        am = a{m+1};
        bm = b{m+1};

        a_next = (A * am - B * bm) / Nm;
        b_next = (A * bm + B * am) / Nm;

        % Re-orthogonalize lightly against previously built columns
        % to control numerical drift
        for k = 1:m
            a_next = a_next - (a{k+1}' * a_next) * a{k+1} ...
                            - (b{k+1}' * a_next) * b{k+1};
            b_next = b_next - (a{k+1}' * b_next) * a{k+1} ...
                            - (b{k+1}' * b_next) * b{k+1};
        end
        a_next = a_next - (a0' * a_next) * a0;
        b_next = b_next - (a0' * b_next) * a0;

        % Normalize so that [sqrt(2)a_m, sqrt(2)b_m] becomes orthonormal
        sa = norm(a_next);
        sb = norm(b_next);

        if sa < 1e-14 || sb < 1e-14
            error('Degenerate recurrence at m = %d. Check generator convention.', m+1);
        end

        % Average the two norms to preserve pair symmetry
        s = 0.5 * (sa + sb);
        a_next = a_next / s;
        b_next = b_next / s;

        % Enforce orthogonality inside the pair
        b_next = b_next - (a_next' * b_next) / (a_next' * a_next) * a_next;
        b_next = b_next / norm(b_next) * norm(a_next);

        if b_next' * Lu * a_next > 0
            b_next = -b_next;
        end

        a{m+2} = a_next;
        b{m+2} = b_next;
    end

    % Assemble F in the order:
    % [sqrt(2)a_l, sqrt(2)b_l, ..., sqrt(2)a_1, sqrt(2)b_1, a_0]
    F = zeros(n, n);
    col = 1;

    %for m = l:-1:1
    %    F(:, col) = sqrt(2) * a{m+1};
    %    col = col + 1;
    %    F(:, col) = sqrt(2) * b{m+1};
    %    col = col + 1;
    %end
    %
    %F(:, col) = a0;

    F(:,l+1) = a0;
    for m = l:-1:1
        F(:,l+1-m) = sqrt(2) * a{m+1};
        F(:,l+1+m) = sqrt(2) * b{m+1};
    end

    % Final orthogonalization/polar correction for numerical robustness
    [Q, ~] = qr(F, 0);
    F = Q;

    % Optional sign fix: enforce last column to align with a0
    if F(:, end)' * a0 < 0
        F(:, end) = -F(:, end);
    end
end