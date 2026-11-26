function [Out] = rsh_rot_ivanic_u(m, n, l, R, Mprev)
%RSH_ROT_IVANIC_U  Function U of the Ivanic-Ruedenberg recursion.
%
%   Out = RSH_ROT_IVANIC_U(m, n, l, R, Mprev)
%
%   Helper of RSH_ROT_IVANIC (Ivanic & Ruedenberg, J. Phys. Chem. 1996/1998).
%   l is the band, R the 3x3 band-1 matrix, Mprev the matrix of band l-1.
    Out = rsh_rot_ivanic_p(0, m, n, l, R, Mprev); 
end
