function [Ml, M] = rsh_rot_lie_eig(lmax, x, gens)
%RSH_ROT_LIE_EIG  RSH rotation matrices of bands 0..lmax by anti-diagonal factorization.
%
%   [Ml, M] = rsh_rot_lie_eig(lmax, x)
%   [Ml, M] = rsh_rot_lie_eig(lmax, x, gens)
%
%   x    : Rodrigues vector (theta * unit axis)
%   gens : optional cell array, gens{l} = struct('Lx',..,'Ly',..,'Lz',..) for l = 1..lmax.
%          Pass it to avoid rebuilding the generators at every call.
%   Ml   : (2lmax+1)x(2lmax+1) matrix of band lmax
%   M    : cell array, M{l+1} = M_l, l = 0..lmax
%
%   M_l = F_l * M_lz(theta) * F_l', with F_l from rsh_rot_to_z.
%
%   This is the method proposed in the paper (called "Lie-AD" in the
%   experiments): M_l(theta, u) = expm(theta * u.L) is evaluated without any
%   matrix exponential or complex eigen-decomposition. F_l is a real
%   orthogonal matrix that depends only on the axis u, and M_lz(theta) is the
%   elementary rotation about z, whose only non-zero entries are cos(m theta)
%   on the diagonal and +/- sin(m theta) on the anti-diagonal.
%
%   Example
%       x = 0.8 * [1; 2; 2] / 3;              % 0.8 rad about [1 2 2]/3
%       [M8, M] = rsh_rot_lie_eig(8, x);      % bands 0..8
%       [Lx, Ly, Lz] = rsh_lie_generator(8);
%       norm(M8 - expm(x(1)*Lx + x(2)*Ly + x(3)*Lz))   % ~1e-15
%
%   See also RSH_ROT_TO_Z, RSH_LIE_GENERATOR, RSH_ROT_LIE, RSH_ROT_COEFFS.
    x = x(:);
    theta = norm(x);
    M = cell(1, lmax+1); M{1} = 1; Ml = 1;
    if theta < 1e-14                                % identity rotation
        for l = 1:lmax, M{l+1} = eye(2*l+1); end
        Ml = M{end};
        return;
    end
    x_hat = x / theta;
    Mz = 1; cl = 1; sl = 0;
    ct = cos(theta); st = sin(theta);
    for l = 1:lmax
        % cos(l theta), sin(l theta) by angle addition
        cl2 = cl * ct - sl * st;
        sl  = sl * ct + cl * st;
        cl  = cl2;
        vz = zeros(2*l-1, 1);
        Mz = [cl, vz', sl;
              vz, Mz,  vz;
             -sl, vz', cl];
        if nargin >= 3 && ~isempty(gens)
            Fl = rsh_rot_to_z(l, x_hat, gens{l});
        else
            Fl = rsh_rot_to_z(l, x_hat);
        end
        Ml = Fl * Mz * Fl';
        M{l+1} = Ml;
    end
end
