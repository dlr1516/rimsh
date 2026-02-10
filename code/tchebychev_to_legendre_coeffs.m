% [Coeffs] = tchebychev_poly_coeffs(n)
% Returns a matrix with the coefficients of Tchebychev polynomials.
% The rows from 1 to (n+1) 
function [Coeffs2] = tchebychev_to_legendre_coeffs(n)
  if (n==0)
      Coeffs = eye(1, floor(n/2)+1);
      Coeffs2 = eye(1, floor(n/2)+1);
      return
  end

  Coeffs = zeros(1, floor(n/2)+1);
  for k=0:floor(n/2)
    Coeffs(k+1) = factor_up(-0.5,k) * n * factor_up(1,n-k-1) * (n-2*k+0.5) / ...
                  (2 * factor_up(1,k) * factor_up(0.5,n-k+1));
  end

  % Alternative computation avoiding large number at numerator and denominator
  Coeffs2 = zeros(1,floor(n/2)+1);
  % Computes the initial item of the sequence
  Coeffs2(1) = (n+0.5);
  for j=1:n
    Coeffs2(1) = Coeffs2(1) * (2 * j) / (2 * j + 1);
  end
  for k=1:floor(n/2)
    Coeffs2(k+1) = Coeffs2(k) * ...
       (2*k-3) * (n-2*k+0.5) * (n-k+1.5) / ( (2*k) * (n-2*k+2.5) * (n-k) );
  end
  Comparison = [Coeffs', Coeffs2'];
end 
