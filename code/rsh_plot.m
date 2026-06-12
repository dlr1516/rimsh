function [Z, latDeg, lonDeg,Theta,Phi] = rsh_plot(rsh_coeffs,res,title_name,grid_mode)
  if (~exist("res"))
      res = 180;
  end
  if (~exist("title_name"))
      grid_mode = "no title";
  end
  if (~exist("grid_mode"))
      grid_mode = "none";
  end

  lmax = floor(sqrt(length(rsh_coeffs))) - 1;
  if ((lmax+1)^2 ~= length(rsh_coeffs))
      disp(['Invalid size of SH coefficients vector ' num2str(length(rsh_coeffs)) ...
          ': it must be the square of interger like ' num2str((lmax+1)^2)]);
      error((lmax+1)^2 ~= length(rsh_coeffs));
  end
  step = pi / res;
  theta = 0:step:pi;
  phi = 0:step:2*pi;
  [Yreal,Theta,Phi] = rsh_eval(lmax,theta,phi);
  latDeg = 90 - 180/pi*theta;
  lonDeg = 180/pi*phi;
  % 1) Yreal size is [(lmax+1)^2, thetaNum, phiNum]
  % 2) rsh_coeffs size is (lmax+1)^2
  % 3) Z is the tensor product on contraction dimension 1 of Yreal and 
  %    1 of rsh_coeffs
  %    -> Z(theta,phi) has size [thetaNum, phiNum]
  Z = tensorprod(Yreal, rsh_coeffs,1,1);
  %size_Rsh = size(Z)
  pc = pcolor(lonDeg,latDeg,Z);
  if grid_mode == "none"
      set(pc, 'EdgeColor', 'none');
  end
  title(title_name)
  xlabel('lon [deg]')
  ylabel('lat [deg]')
  colorbar()
end