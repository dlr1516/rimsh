% [Ml, M] = rsh_rot_lie_eig(lmax, x)
%
% Computes the rotation matrices of bands l=0, ..., lmax of real 
% spherical harmonics (RSH) corresponding to the Lie algebra so(3) 
% paramaters x.
%
% Input: 
% - lmax: maximum order l of RSH (integer);
% - x: 3-vector of Lie algrebra parameters.
%
% Output:
% - Ml: the RSH rotation matrix for band lmax with size 2*lmax+1
% - M: cell array of all the RSH rotation matrices for bands from 0 to
%      lmax and the last one M{lmax} is a copy of M. 
function [Ml, M, L] = rsh_rot_lie_eig(lmax, x)
   x = x(:);
   theta = norm(x);
   if theta < 1e-14
       x_hat = zeros(3,1);
   else
       x_hat = x ./ theta;
   end

   M = {1};
   Mz = 1;
   ct = cos(theta);
   st = sin(theta);
   cl = 1;
   sl = 0;
   for l=1:lmax
     cl2 = cl * ct - sl * st;
     sl2 = sl * ct + cl * st;
     cl = cl2;
     sl = sl2;
     vz = zeros(2*l-1,1);
     Mz = [cl, vz', sl;
           vz, Mz,  vz;
          -sl, vz', cl];
    Fl = rsh_rot_to_z(l, x_hat);
    Ml = Fl * Mz * Fl';
    M{l+1} = Ml;
  end
end