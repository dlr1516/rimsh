function [Yrsh,Theta,Phi] = rsh_eval(lmax, theta, phi)
% RSH_EVAL
% It computes the value of Real Spherical Harmonics (RSH). 
% It can be called in two ways according to the number of arguments.
%
% [Yrsh, Theta, Phi] = rsh_eval(lmax, theta, phi)
% Computes the Real Spherical Harmonics (RSH) for indices l=0, 1, ..., lmax 
% at the given values of theta and phi where theta and phi can be column vectors.
%
% INPUTS:
% - lmax: maximum band size of RSH
% - theta: vector of colatitude coordinates of size nt x 1
%     nt = length(theta)
% - phi: vector of azimuth coordinates of size np x 1
%     np = length(phi)
%
% OUTPUTS:
% - Yrsh: tensor of RSH values with concatenation of all the bands 
%     l=0, 1, ..., lmax of size (lmax + 1)^2 x nt x np
% - Theta: matrix with repeated value of theta in each row
%     with matrix size nt x np
% - Phi: matrix with repeated value of phi in each column
%     with matrix size nt x np
%
% [Yrsh, Theta, Phi] = rsh_eval(lmax, vec) 
% Computes the Real Spherical Harmonics (RSH) for indices l=0, 1, ..., lmax 
% at the given unitary vector vec using its polar coordinates:
%   theta = atan2(sqrt(vec(1,:).^2 * vec(2,:).^2), vec(3,:))
%   phi = atan2(v)
%
% INPUTS:
% - lmax: maximum band size of RSH
% - vec: matrix of size 3 x nv where each column is a unit vector 
%
% OUTPUTS:
% - Yrsh: tensor of RSH values with concatenation of all the bands 
%     l=0, 1, ..., lmax of size (lmax + 1)^2 x nv x nv
% - Theta: matrix with repeated value of theta in each row
%     with matrix size nv x nv
% - Phi: matrix with repeated value of phi in each column
%     with matrix size nv x nv
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
  Yrsh=zeros(sh_lm_to_index(lmax,lmax),thetaNum,phiNum);
  for l=0:lmax
    CosSin = [sin((l:-1:1)' * phi); ones(1,phiNum); cos((1:l)' * phi)];
    CosSin = reshape(CosSin, [(2*l+1) 1 phiNum]);
    CosSin = repmat(CosSin, [1 thetaNum 1]);

    % legendre(l,cos(theta),'sch') returns a vector of size [(l+1) length(theta)]
    Pl=sqrt((2*l+1)/(4*pi))*legendre(l,cos(theta),'sch');
    Ptot = [Pl(l+1:-1:2,:); Pl(1,:); Pl(2:l+1,:)];
    Ptot = reshape(Ptot, [(2*l+1), thetaNum, 1]);
    Ptot = repmat(Ptot, [1 1 phiNum]);
    %size_Ptot = size(Ptot)

    Yrsh(sh_lm_to_index(l,-l:l),:,:) = Ptot .* CosSin;
  end
end
