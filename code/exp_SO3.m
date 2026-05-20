function [R] = exp_SO3(x)
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