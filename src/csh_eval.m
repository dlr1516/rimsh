function [Ycsh,Theta,Phi] = csh_eval(lmax, theta, phi)
%CSH_EVAL  Evaluate complex spherical harmonics of bands 0..lmax.
%
% It computes the value of Complex Spherical Harmonics (CSH). 
% It can be called in two ways according to the number of arguments.
%
% [Ycsh, Theta, Phi] = csh_eval(lmax, theta, phi)
% Computes the Complex Spherical Harmonics (CSH) for indices l=0, 1, ..., lmax 
% at the given values of theta and phi where theta and phi can be column vectors.
%
% INPUTS:
% - lmax: maximum band (degree) l
% - theta: vector of colatitude coordinates of size nt x 1
%     nt = length(theta)
% - phi: vector of azimuth coordinates of size np x 1
%     np = length(phi)
%
% OUTPUTS:
% - Ycsh: tensor of CSH values with concatenation of all the bands 
%     l=0, 1, ..., lmax of size (lmax + 1)^2 x nt x np
% - Theta: matrix with repeated value of theta in each row
%     with matrix size nt x np
% - Phi: matrix with repeated value of phi in each column
%     with matrix size nt x np
%
% [Ycsh, Theta, Phi] = csh_eval(lmax, vec) 
% Computes the Complex Spherical Harmonics (CSH) for indices l=0, 1, ..., lmax 
% at the given unitary vector vec using its polar coordinates:
%   theta = atan2(sqrt(vec(1,:).^2 + vec(2,:).^2), vec(3,:))
%   phi = atan2(vec(2,:), vec(1,:))
%
% INPUTS:
% - lmax: maximum band (degree) l
% - vec: matrix of size 3 x nv where each column is a unit vector 
%
% OUTPUTS:
% - Ycsh: tensor of CSH values with concatenation of all the bands 
%     l=0, 1, ..., lmax of size (lmax + 1)^2 x nv x nv
% - Theta: matrix with repeated value of theta in each row
%     with matrix size nv x nv
% - Phi: matrix with repeated value of phi in each column
%     with matrix size nv x nv
% NOTE: theta and phi of all the nv vectors are combined on a grid, so the
% values at the nv input points are Y(:, k, k). For a single vector (nv = 1)
% Y is simply a (lmax+1)^2 x 1 vector.
%
% Convention: orthonormal on the unit sphere, no Condon-Shortley phase,
%   Y_l^m = sqrt((2l+1)/(4pi) (l-|m|)!/(l+|m|)!) P_l^|m|(cos theta) exp(i m phi),
% hence Y_l^{-m} = conj(Y_l^m).
%
% Coefficients are packed band after band, m = -l..l (see SH_LM_TO_INDEX).
%
% DEPENDENCIES:
% This function uses Matlab function legendre() for computation 
% of associated Legendre polynomials with Schmidt semi-normalization. 
%
  if (nargin == 2)
    vec = theta;
    [dim, num] = size(vec);
    assert(dim == 3);
    theta =  atan2(sqrt(vec(1,:).^2 + vec(2,:).^2), vec(3,:));
    phi = atan2(vec(2,:), vec(1,:));
  end
  thetaNum = length(theta);
  phiNum = length(phi);
  theta= reshape(theta, [thetaNum 1]);
  phi = reshape(phi, [1 phiNum]);
  [Phi,Theta] = meshgrid(phi, theta);
  Ycsh=zeros(sh_lm_to_index(lmax,lmax),thetaNum,phiNum);
  for l=0:lmax
    ExpPhi = exp(1i * (-l:l)' * phi);
    ExpPhi = reshape(ExpPhi, [(2*l+1) 1 phiNum]);
    ExpPhi = repmat(ExpPhi, [1 thetaNum 1]);

    % legendre(l,cos(theta),'sch') returns a vector of size [(l+1) length(theta)]
    % Schmidt semi-normalized P_l^m carry an extra sqrt(2) for m > 0 that is
    % needed by real SH but not by orthonormal complex SH: remove it.
    Pl=sqrt((2*l+1)/(4*pi))*legendre(l,cos(theta),'sch');
    Pl(2:end,:) = Pl(2:end,:) / sqrt(2);
    Ptot = [Pl(l+1:-1:2,:); Pl(1,:); Pl(2:l+1,:)];
    Ptot = reshape(Ptot, [(2*l+1), thetaNum, 1]);
    Ptot = repmat(Ptot, [1 1 phiNum]);

    Ycsh(sh_lm_to_index(l,-l:l),:,:) = Ptot .* ExpPhi;
  end
end
