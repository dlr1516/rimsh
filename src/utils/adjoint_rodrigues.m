function [Omega,A,B] = adjoint_rodrigues(a, b, l)
%ADJOINT_RODRIGUES  Adjoint action exp(A) * B * exp(-A) in closed form.
%
%   [Omega, A, B] = ADJOINT_RODRIGUES(a, b, l)
%
%   A = a(1)*Lx + a(2)*Ly + a(3)*Lz and B = b(1)*Lx + b(2)*Ly + b(3)*Lz are
%   elements of the band-l representation of so(3) (see RSH_LIE_GENERATOR).
%   Since ad_A has eigenvalues 0, +/- i*norm(a), the adjoint map is
%       Omega = B + sin(k)/k [A,B] + (1-cos(k))/k^2 [A,[A,B]],  k = norm(a).
%
%   See also RSH_LIE_GENERATOR.
  [Lx, Ly, Lz] = rsh_lie_generator(l);
  A = Lx * a(1) + Ly * a(2) + Lz * a(3); 
  B = Lx * b(1) + Ly * b(2) + Lz * b(3); 
  ka = norm(a);
  lb_AB = (A * B - B * A);
  lb_AAB = (A * lb_AB - lb_AB * A);
  if ka < 1e-12
    Omega = B + lb_AB;                    % first-order limit
  else
    Omega = B + sin(ka) / ka * lb_AB + (1 - cos(ka)) / ka^2 * lb_AAB;
  end
end
