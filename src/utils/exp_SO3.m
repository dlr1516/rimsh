function [R] = exp_SO3(x)
%EXP_SO3  Exponential map from R^3 (rotation vector) to SO(3).
%
%   R = EXP_SO3(x)
%
%   Rodrigues formula: R = I + sin(t) K + (1 - cos(t)) K^2 with t = norm(x)
%   and K the skew-symmetric matrix of the unit axis x/t.
%
%   See also LOG_SO3.
  assert(length(x) == 3)
  theta = norm(x);
  if theta < 1e-6
    R = eye(3);
  else
    % Unit vector axis
    k = x/theta; 
    K = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
    % Rodrigues formula
    R = eye(3) + sin(theta)*K + (1-cos(theta))*(K^2);
  end
end
