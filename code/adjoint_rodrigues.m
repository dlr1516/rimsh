% [X] = adjoint_rodrigues(a, b)
%
% Computes the adjoint map: Omega = exp(A) * B * exp(-A).
% The matrices A, B are representation of the Lie algebra LA) representation
% of so(3) in (2*l+1) spherical harmonics space. 
% A and B are defined as linear combination of the generators of the LA:
%   A = Lx * a(1) + Ly * a(2) + Lz * a(3)
%   B = Lx * b(1) + Ly * b(2) + Lz * b(3);
% where a and b are 3-vectors. 
function [Omega,A,B] = adjoint_rodrigues(a, b, l)
  %
  [Lx, Ly, Lz] = rsh_lie_generator(l);
  A = Lx * a(1) + Ly * a(2) + Lz * a(3); 
  B = Lx * b(1) + Ly * b(2) + Lz * b(3); 
  ka = norm(a);
  kb = norm(b);
  lb_AB = (A * B - B * A);
  lb_AAB = (A * lb_AB - lb_AB * A);
  Omega = B + sin(ka) /ka * lb_AB + (1-cos(ka)) / (ka^2) * lb_AAB;
end