% [Ml, M] = rsh_rot_lie(lmax, x)
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
function [Ml, M, L] = rsh_rot_lie(lmax, x, data)
  % Creates the Lie algebra generators
  if nargin < 3
    data = struct();
  end
  if ~isfield(data, 'L')
    for l=1:lmax
      Ll = zeros(2*l+1,2*l+1,3);
      [Ll(:,:,1), Ll(:,:,2), Ll(:,:,3)] = rsh_lie_generator(l); 
      data.L{l} = Ll;
    end
  end

  M = {1};
  for l=1:lmax
    [Lx, Ly, Lz] = rsh_lie_generator(l);
    Ml = expm(Lx * x(1) + Ly * x(2) + Lz * x(3));
    %Ml = expm(data.L{1}(:,:,1) * x(1) + expm(data.L{1}(:,:,2) * x(2) + expm(data.L{1}(:,:,3) * x(3));
    M{l+1} = Ml;
  end
end