function [f, g, H] = rsh_corr_grad_hess(lmax, rsh_src, rsh_dst, x)
  [r,c] = size(x);
  if r == 3 & c == 3
    R = x;
    x = log_SO3(R);
  elseif r * c == 3
    R = exp_SO3(x);  
  else
    error('invalid input x: it is either a 3x3 rotation matrix or a 3x1 vector of Lie algebra parameters')
  end

  f = 0;
  g = zeros(3, 1); 
  H = zeros(3, 3); 

  for l=1:lmax
    L = zeros(2*l+1,2*l+1,3);
    [L(:,:,1), L(:,:,2), L(:,:,3)] = rsh_lie_generator(l); 
    % Takes the order l slices of RSH source and destination coefficients
    s = rsh_src(sh_lm_to_index(l,-l:l));
    d = rsh_dst(sh_lm_to_index(l,-l:l));
    M = expm(x(1) * L(:,:,1) + x(2) * L(:,:,2) + x(3) * L(:,:,3));
    dt = M' * d;
    %st = M * s;
    f = f + d' * M * s;
    %disp('[s | dt]')
    %disp([s, dt])
    %disp(['dt^T * Lx * s ' num2str(dt' * L(:,:,1) * s)])
    % Gradient
    for j = 1:3
      g(j) = g(j) + dt' * L(:,:,j) * s;
      %disp(['g(' num2str(j) ') l ' num2str(l) ': incr ' num2str(d' * L(:,:,j) * s)])
    end
    for j = 1:3
      for k = j+1:3
        Inc = dt' * (L(:,:,j) * L(:,:,k) + L(:,:,k) * L(:,:,j)) * s ./ 2;
        H(j,k) = H(j,k) + Inc;
        if j ~= k
           H(k,j) = H(k,j) + Inc;
        end
      end
    end
  end
end