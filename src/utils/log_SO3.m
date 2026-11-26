function w = log_SO3(R)
%LOG_SO3  Logarithm map from SO(3) to R^3 (rotation vector).
%   w = LOG_SO3(R) returns the 3x1 Lie algebra vector w such that
%   R = expm(hat(w)) approximately.
%
%   Input:
%       R : 3x3 rotation matrix
%
%   Output:
%       w : 3x1 rotation vector (axis * angle)
%
%   See also EXP_SO3.

    % Optional: project to nearest rotation matrix for robustness
    [U, ~, V] = svd(R);
    R = U * diag([1 1 det(U*V')]) * V';

    tr = trace(R);
    c = (tr - 1) / 2;
    c = min(1, max(-1, c));   % clamp for numerical safety
    theta = acos(c);

    % Small-angle case
    if theta < 1e-8
        % Since R-R' ~ 2*hat(w)
        w = vee(0.5 * (R - R'));
        return;
    end

    % Near-pi case
    if abs(pi - theta) < 1e-6
        % Use diagonal entries to recover axis
        A = (R + eye(3)) / 2;

        axis = zeros(3,1);
        axis(1) = sqrt(max(A(1,1), 0));
        axis(2) = sqrt(max(A(2,2), 0));
        axis(3) = sqrt(max(A(3,3), 0));

        % Fix signs using off-diagonal entries
        if A(1,2) < 0, axis(2) = -axis(2); end
        if A(1,3) < 0, axis(3) = -axis(3); end

        % If axis is nearly zero due to numerical issues, use fallback
        if norm(axis) < 1e-8
            [~, idx] = max(diag(A));
            axis = A(:, idx);
            axis = axis / norm(axis);
        else
            axis = axis / norm(axis);
        end

        w = theta * axis;
        return;
    end

    % Generic case
    W = (theta / (2*sin(theta))) * (R - R');
    w = vee(W);
end

function W = hat(w)
    W = [   0   -w(3)  w(2);
          w(3)   0    -w(1);
         -w(2)  w(1)   0   ];
end

function w = vee(W)
    w = [W(3,2); W(1,3); W(2,1)];
end
