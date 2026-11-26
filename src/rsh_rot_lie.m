function [Ml, M] = rsh_rot_lie(lmax, x, gens)
%RSH_ROT_LIE  RSH rotation matrices of bands 0..lmax by matrix exponential.
%
%   [Ml, M] = RSH_ROT_LIE(lmax, x)
%   [Ml, M] = RSH_ROT_LIE(lmax, x, gens)
%
%   Reference (generic) implementation of the Lie-algebraic formula
%       M_l(x) = expm(x(1)*Lx + x(2)*Ly + x(3)*Lz),
%   with the generators of RSH_LIE_GENERATOR. For the fast closed-form
%   evaluation proposed in the paper use RSH_ROT_LIE_EIG.
%
%   Inputs
%     lmax : maximum band (non-negative integer)
%     x    : 3-vector, Rodrigues / rotation vector (theta * unit axis)
%     gens : optional cell array, gens{l} = struct('Lx',..,'Ly',..,'Lz',..)
%            for l = 1..lmax, to avoid rebuilding the generators at each call
%
%   Outputs
%     Ml   : (2lmax+1)x(2lmax+1) rotation matrix of band lmax
%     M    : cell array with M{l+1} = M_l, l = 0..lmax
%
%   See also RSH_ROT_LIE_EIG, RSH_LIE_GENERATOR, RSH_ROT_COEFFS.

  M = cell(1, lmax+1);
  M{1} = 1;
  Ml = 1;
  for l = 1:lmax
    if nargin >= 3 && ~isempty(gens)
      Lx = gens{l}.Lx; Ly = gens{l}.Ly; Lz = gens{l}.Lz;
    else
      [Lx, Ly, Lz] = rsh_lie_generator(l);
    end
    Ml = expm(Lx * x(1) + Ly * x(2) + Lz * x(3));
    M{l+1} = Ml;
  end
end
