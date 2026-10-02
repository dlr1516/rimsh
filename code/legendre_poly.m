function [y] = legendre_poly(n, x)
  yprev2 = x;
  yprev1 = 1.5 .* x.^2 - 0.5;
  if n == 0
    y = ones(size(x));
  elseif n==1
    y = yprev2;
  elseif n==2
    y = yprev1;            
  else 
    for l=3:n
      y = ((2.*l-1) .* x .* yprev1 - (l-1) .* yprev2) ./ l;
      yprev2 = yprev1;
      yprev1 = y;
    end
  end
return