function [Z, latDeg, lonDeg, Theta, Phi] = rsh_plot(rsh_coeffs, res, title_name, grid_mode)
%RSH_PLOT  Plot a function given by its RSH coefficients on a lat/lon map.
%
%   RSH_PLOT(c)
%   RSH_PLOT(c, res, title_name, grid_mode)
%   [Z, latDeg, lonDeg, Theta, Phi] = RSH_PLOT(...)
%
%   Inputs
%     c          : (lmax+1)^2 vector of RSH coefficients (see SH_LM_TO_INDEX)
%     res        : number of samples on [0, pi] (grid step pi/res), default 180
%     title_name : plot title, default ''
%     grid_mode  : 'none' (default) hides the cell edges of PCOLOR, any other
%                  value keeps them
%
%   Outputs
%     Z          : function values on the (colatitude x azimuth) grid
%     latDeg     : latitude [deg] of the grid rows
%     lonDeg     : longitude [deg] of the grid columns
%     Theta, Phi : grids of colatitude and azimuth [rad]
%
%   See also RSH_EVAL.

  if nargin < 2 || isempty(res),        res = 180; end
  if nargin < 3 || isempty(title_name), title_name = ''; end
  if nargin < 4 || isempty(grid_mode),  grid_mode = 'none'; end

  lmax = round(sqrt(numel(rsh_coeffs))) - 1;
  if (lmax+1)^2 ~= numel(rsh_coeffs)
    error('rsh_plot:size', ...
          'length of rsh_coeffs (%d) must be (lmax+1)^2.', numel(rsh_coeffs));
  end
  step = pi / res;
  theta = 0:step:pi;
  phi = 0:step:2*pi;
  [Yreal, Theta, Phi] = rsh_eval(lmax, theta, phi);
  latDeg = 90 - 180/pi*theta;
  lonDeg = 180/pi*phi;
  % Yreal is (lmax+1)^2 x nTheta x nPhi: contract the first dimension with c
  [~, nt, np] = size(Yreal);
  Z = reshape(rsh_coeffs(:).' * reshape(Yreal, [], nt*np), nt, np);
  pc = pcolor(lonDeg, latDeg, Z);
  if strcmpi(grid_mode, 'none')
    set(pc, 'EdgeColor', 'none');
  end
  title(title_name)
  xlabel('lon [deg]')
  ylabel('lat [deg]')
  colorbar();
end
