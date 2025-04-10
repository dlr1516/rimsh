function [Ml, M] = rsh_rot_ivanic(R, lmax)
% RSH_ROT_IVANIC
% Comutes the rotation matrix of the Real Spherical Harmonics (RSH)
% coeffients according to the formula given in: 
%
% Ivanic, Ruedenberg, "Rotation matrices for real spherical harmonics. 
% Direct determination by recursion", Journal of Physical Chemistry, 
% 1995, 100(15), 6342-6347.
%
% INPUTS:
% - R: the 3D rotation matrix applied to the reference frame of SH
% - lmax: maximum band of spherical harmonics
%
% OUTPUTS:
% - Ml: the RSH rotation matrix for band lmax with size 2*lmax+1
% - M: cell array of all the RSH rotation matrices for bands from 0 to
%      lmax and the last one M{lmax} is a copy of Ml 
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
    u_mat = zeros(2*l+1,2*l+1);
    v_mat = zeros(2*l+1,2*l+1);
    w_mat = zeros(2*l+1,2*l+1);
    for m = -l:l
      for n = -l:l
        [u, v, w] = rsh_rot_ivanic_coeffs(m, n, l);
        u_mat(m+l+1,n+l+1) = u;
        v_mat(m+l+1,n+l+1) = v;
        w_mat(m+l+1,n+l+1) = w;
        %disp(['l ' num2str(l) ' row ' num2str(m) ' col ' num2str(n) ': u ' num2str(u) ' v ' num2str(v) ' w ' num2str(w)])
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
    %startIdx = l^2 + 1;
  end
end
