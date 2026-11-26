function [Ml, M] = rsh_rot_ivanic(R, lmax)
%RSH_ROT_IVANIC  RSH rotation matrices by the Ivanic-Ruedenberg recursion.
%
%   [Ml, M] = RSH_ROT_IVANIC(R, lmax)
%
% Computes the rotation matrix of the Real Spherical Harmonics (RSH)
% coefficients according to the formula given in: 
%
% Ivanic, Ruedenberg, "Rotation matrices for real spherical harmonics. 
% Direct determination by recursion", Journal of Physical Chemistry, 
% 1996, 100(15), 6342-6347.
%
% and their errata corrige:
%
% Ivanic, Ruedenberg, "Rotation matrices for real spherical harmonics. 
% Direct determination by recursion", Journal of Physical Chemistry,
% 1998, 102(45), 9099-9100.
%
% INPUTS:
% - R: the 3D rotation matrix applied to the reference frame of SH
% - lmax: maximum band of spherical harmonics
%
% OUTPUTS:
% - Ml: the RSH rotation matrix for band lmax with size 2*lmax+1
% - M: cell array of all the RSH rotation matrices for bands from 0 to
%      lmax (M{l+1} = M_l); M{end} is a copy of Ml.
%
% NOTE: Reference baseline of the paper. With the conventions of this library
% rsh_rot_ivanic(R, lmax) gives the same matrices as
% rsh_rot_lie_eig(lmax, log_SO3(R)) (checked in tests/run_tests.m).
%
% See also RSH_ROT_LIE_EIG, RSH_ROT_LIE.
  EPS = 1e-6;
  % M1 = [R(2,2), R(2,3), R(2,1);  
  %     R(3,2), R(3,3), R(3,1);
  %     R(1,2), R(1,3), R(1,1)];
  E = [0 0 1;
       1 0 0;
       0 1 0];
  Ml = E' * R * E;
  M = {1, Ml};
  for l=2:lmax
    % This function also works with symbolic matrix R   
    if strcmpi(class(R),'sym')
      Ml = sym(zeros(2*l+1,2*l+1));
    else
      Ml = zeros(2*l+1,2*l+1);
    end
    for m = -l:l
      for n = -l:l
        [u, v, w] = rsh_rot_ivanic_coeffs(m, n, l);
        Ml(m+l+1,n+l+1) = 0.0;
        if (abs(u) > EPS)
          Ml(m+l+1,n+l+1) = Ml(m+l+1,n+l+1) + u * rsh_rot_ivanic_u(m, n, l, M{2}, M{end});
        end
        if (abs(v) > EPS)
          Ml(m+l+1,n+l+1) = Ml(m+l+1,n+l+1) + v * rsh_rot_ivanic_v(m, n, l, M{2}, M{end});
        end
        if (abs(w) > EPS)
          Ml(m+l+1,n+l+1) = Ml(m+l+1,n+l+1) + w * rsh_rot_ivanic_w(m, n, l, M{2}, M{end});
        end
      end
    end
    M{end+1} = Ml;
  end
end
