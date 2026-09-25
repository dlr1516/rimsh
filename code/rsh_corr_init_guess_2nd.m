function [x_init, R_init] = rsh_corr_init_guess_2nd(lmax, rsh_src, rsh_dst)
% RSH_CORR_INIT_GUESS
% It estimates the value of the rotation that better aligns two functions
% on sphere expanded in Real Spherical Harmonics (RSH). 
% The estimates is meant as an initial guess for a local search optimization
% such as Gauss-Newton method. 
%
% MOTIVATION
% .....
%
% SUGGESTION
% The method operates on the coeffients of a single band l. 
% The user must select  
%
% INPUTS:
% - l: the band size of RSH to be used in optmization
% - rsh_src: RSH coefficients of source function
% - rsh_dst: RSH coefficients of destination function
%
% OUTPUTS:
% - x_init: Rodrigues paramaters of rotation, i.e. the axis vector multiplied 
%           by the rotation angle in radians;
% - R_init: the the initial guess in the form of rotation matrix. 
%
  Hs = zeros(3,3);
  Hd = zeros(3,3);
  RotAxes = zeros(3,0);
  for l=1:lmax    % TODO: start with single band l do the same for l=0:lmax and sum teh result
    L = zeros(2*l+1,2*l+1,3);
    [L(:,:,1), L(:,:,2), L(:,:,3)] = rsh_lie_generator(l);
    s = rsh_src(sh_lm_to_index(l,-l:l));
    d = rsh_dst(sh_lm_to_index(l,-l:l));
    IncHs = zeros(3,3);
    IncHd = zeros(3,3);
    for ii=1:3
      for jj=ii:3
        IncHs(ii,jj) = s' * (L(:,:,ii) * L(:,:,jj) + L(:,:,jj) * L(:,:,ii)) * s;
        IncHs(jj,ii) = IncHs(ii,jj);
        IncHd(ii,jj) = d' * (L(:,:,ii) * L(:,:,jj) + L(:,:,jj) * L(:,:,ii)) * d;
        IncHd(jj,ii) = IncHd(ii,jj);
      end
    end
    Hs = Hs + IncHs;
    Hd = Hd + IncHd;
    if mod(l,2) == 0
      a = [d' * L(:,:,1) * s; d' * L(:,:,2) * s; d' * L(:,:,3) * s];
      a = a / norm(a);
      RotAxes = [RotAxes, a];
    end
  end
  [Vs, Ds] = eig(Hs);
  [Vd, Dd] = eig(Hd);
  disp('diag(Ds): ')
  disp(diag(Ds)')
  disp('diag(Dd): ')
  disp(diag(Dd)')

  R_init = Vs' * Vd;   % possible reflection to be removed: check eigenvectors
  x_init = log_SO3(R_init);

  disp('R_init')
  disp(R_init)
  disp(['x_init: angle[deg] ' num2str(180/pi * norm(x_init)) ' axis'])
  disp(x_init / norm(x_init)) 

  disp('RotAxes')
  disp(RotAxes)
return
