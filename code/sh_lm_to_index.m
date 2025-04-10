function [idx] = sh_lm_to_index(l,m)
  if (abs(m) > l) 
      error(['invalid index m ' num2str(m) ': it must be between -l ' num2str(-l) ' and l ' num2str(l)])
  end
  idx = l * l + m + l + 1;
end
