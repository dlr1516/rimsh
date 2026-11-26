function [Out] = rsh_rot_ivanic_p(i, a, b, l, R, Mprev)
%RSH_ROT_IVANIC_P  Function P of the Ivanic-Ruedenberg recursion.
%
%   Out = RSH_ROT_IVANIC_P(i, a, b, l, R, Mprev)
%
%   Helper of RSH_ROT_IVANIC (Ivanic & Ruedenberg, J. Phys. Chem. 1996/1998).
%   l is the band, R the 3x3 band-1 matrix, Mprev the matrix of band l-1.
  if (b == l)
    Out = sh_rot_index(R, i, 1) * ...
          sh_rot_index(Mprev, a, l - 1) - ...
          sh_rot_index(R, i, -1) * ...
          sh_rot_index(Mprev, a, -l + 1);
  elseif (b == -l)
    Out = sh_rot_index(R, i, 1) * ...
          sh_rot_index(Mprev, a, -l + 1) + ...
          sh_rot_index(R, i, -1) * ...
          sh_rot_index(Mprev, a, l - 1);
  else
    Out = sh_rot_index(R, i, 0) * sh_rot_index(Mprev, a, b);
  end
end
