function resetRands(seed)
%RESETRANDS  Reset the random number generator to SEED (default 0).
%
%   resetRands
%   resetRands(seed)
%
%   Seeds the generator shared by rand, randn, randi and randperm through
%   RNG (Mersenne Twister), so that the following random draws are
%   reproducible. Any real SEED is accepted: it is converted to the integer
%   range required by RNG, mod(floor(abs(seed)), 2^32), so that, e.g., a
%   time stamp such as posixtime(datetime('now')) can be passed directly.
%
%   NOTE: before this version the function used the legacy calls
%   rand('state', seed) and randn('state', seed), which select MATLAB's old
%   generators. The sequences obtained with the same seed are therefore
%   different from those of the old version.
%
%   Origin: adapted from resetRands.m of the codemeta_utils MATLAB library,
%   as distributed with TPGO-RS (codemeta_utils/math/utilities/resetRands.m),
%   https://github.com/ErnestF22/TPGO-RS
%   The original header carried an %%AUTORIGHTS%% license placeholder.
%   Rewritten here with the RNG interface.
%
%   See also RNG.

  if nargin < 1 || isempty(seed)
    seed = 0;
  end
  if ~isnumeric(seed) || ~isscalar(seed) || ~isfinite(seed) || ~isreal(seed)
    error('resetRands:seed', 'seed must be a finite real scalar.');
  end

  rng(mod(floor(abs(double(seed))), 2^32), 'twister');
end
