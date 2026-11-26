function [idx] = sh_lm_to_index(l,m)
%SH_LM_TO_INDEX  Linear index of the spherical harmonic (l, m) in a packed vector.
%
%   idx = SH_LM_TO_INDEX(l, m)
%
%   Coefficient vectors of bands 0..lmax are packed band after band, each band
%   ordered as m = -l, ..., l, so that idx = l^2 + l + m + 1 (1-based).
%   m can be a vector, e.g. SH_LM_TO_INDEX(l, -l:l) returns the indices of
%   band l. The total length of a vector with bands 0..lmax is (lmax+1)^2.
  if (abs(m) > l) 
      error(['invalid index m ' num2str(m) ': it must be between -l ' num2str(-l) ' and l ' num2str(l)])
  end
  idx = l * l + m + l + 1;
end
