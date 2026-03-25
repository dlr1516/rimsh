function [g,H] = rsh_corr_grad_hess(rsh_src, rsh_dst, lmax)
  g = zeros(3, 1); 
  H = zeros(3, 3); 

  for l=1:lmax
    L = zeros(2*l+1,2*l+1,3);
    [L(:,:,1), L(:,:,2), L(:,:,3)] = rsh_lie_generator(l); 
    % Takes the order l slices of RSH source and destination coefficients
    s = rsh_src(sh_lm_to_index(l,-l:l));
    d = rsh_dst(sh_lm_to_index(l,-l:l));
    %disp('[src | dst]')
    %disp([s, d])
    % Gradient
    for j = 1:3
      g(j) = g(j) + d' * L(:,:,j) * s;
      %disp(['g(' num2str(j) ') l ' num2str(l) ': incr ' num2str(d' * L(:,:,j) * s)])
    end
    for j = 1:3
      for k = j+1:3
        Inc = d' * (L(:,:,j) * L(:,:,k) * L(:,:,k) * L(:,:,j)) * s ./ 2;
        H(j,k) = H(j,k) + Inc;
        if j ~= k
           H(k,j) = H(k,j) + Inc;
        end
      end
    end
  end
end