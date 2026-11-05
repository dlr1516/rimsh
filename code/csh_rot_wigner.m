function [Ml, M] = csh_rot_wigner(x, lmax)
%CSH_ROT_WIGNER  Rotation matrix (Wigner D) for complex spherical harmonics.
%
%   [Ml, M] = CSH_ROT_WIGNER(x, lmax)
%
%   x    : rotation, either
%            - a 3x3 rotation matrix R (orthogonal, det = +1), or
%            - a 3-vector of ZYZ Euler angles (alpha, beta, gamma) in radians,
%              meaning R = Rz(alpha) * Ry(beta) * Rz(gamma).
%   lmax : maximum band index (non-negative integer scalar).
%
%   Ml   : (2*lmax+1) x (2*lmax+1) square rotation matrix acting on the
%          coefficients of the complex spherical harmonics.  Band l occupies
%          rows and columns  lmax-l+1 : lmax+l+1, and inside that block
%          row index = new m-prime, column index = old m, with
%              m-prime, m = -l, -l+1, ..., 0, ..., l
%          (so vector index  m+l+1  for quantum number m).
%
%              Y_l^{m'}(R^{-1} r) = sum_{m=-l}^{l} Ml(m'+l+1, m+l+1) Y_l^{m}(r)
%              f_l^{m'}           = sum_{m=-l}^{l} Ml(m'+l+1, m+l+1) f_l^{m}
%
%   M    : 1 x (lmax+1) cell array with M{l+1} = D^l, the (2*l+1) x (2*l+1)
%          Wigner D matrix of band l.  Hence Ml = M{lmax+1}, but note that
%          the band index m runs over 1:2l+1 = m+l+1 inside M{l}.
%
%   Convention
%   ----------
%   * Complex spherical harmonics with Condon-Shortley phase:
%         Y_l^m(theta,phi) = (-1)^m sqrt((2l+1)/(4 pi) * (l-m)!/(l+m)!)
%                            * P_l^m(cos theta) * exp(i m phi),
%     with P_l^m the associated Legendre function including the (-1)^m factor
%     (that is MATLAB's  legendre(l, cos(theta))  with 'unnorm').
%   * D^l_{m'm}(alpha,beta,gamma) = exp(-i m' alpha) d^l_{m'm}(beta) exp(-i m gamma)
%     (Jacobi convention, the one that is an SO(3) representation matrix).
%   * Active (alibi) rotation: the function / pattern is rotated by R.
%
%   Properties (all hold to round-off)
%   ----------------------------------
%       Ml * Ml'  = I                       (unitary / orthogonal)
%       Ml(x^-1)  = Ml(x)^-1 = Ml(x)'
%       Ml(R1*R2) = Ml(R1) * Ml(R2)         (representation property)
%       Ml is block diagonal, one block per l
%
%   Physics convention: the "physics" D (Akhiezer / Varshalovich, used in much
%   of the angular-momentum literature) is the complex conjugate of the above,
%   D^l_{m'm}^{phys} = conj(D^l_{m'm}), i.e.  M_phys{l} = conj(M{l});  the two
%   differ only in the sign of the elements with odd m', and never for l = 0.
%   If your Y_l^m carries no (-1)^m phase, apply that conjugation.
%
%   Method
%   ------
%   The reduced matrix d^l(beta) is evaluated from its finite Jacobi-polynomial
%   series, which is exact and needs no toolbox:
%
%       d^l_{m'm}(beta) = (-1)^(m-m') sqrt( (l+m)!(l-m)! / ((l+m')!(l-m')!) )
%                         * (cos(beta/2))^(m+m') (sin(beta/2))^(m'-m)
%                         * P_{l-m}^{(m-m', m+m')}(cos(beta)),
%
%   with the Jacobi polynomial expanded in its (l-m+1)-term series.  This is
%   the same expression Varshalovich et al., "Quantum Theory of Angular
%   Momentum", Sec. 4.3, give for d^l.  Cost is O(l^3) per band and the result
%   is exact to round-off for l up to a few tens.
%
%   Example
%   -------
%       lmax = 5;
%       [Ml, M] = csh_rot_wigner([0.3 1.1 -0.7], lmax);   % ZYZ angles
%       l  = 2;
%       f  = randn(2*l+1,1);          % coefficients, ordered m = -l ... l
%       fr = M{l+1} * f;              % rotated coefficients
%       % or from the big matrix:
%       fr = Ml(lmax-l+1:lmax+l+1, lmax-l+1:lmax+l+1) * f;
%
%   See also CSH_ROT_WIGNER_DEMO.

%% ---- argument checks -------------------------------------------------
if ~isnumeric(lmax) || ~isscalar(lmax) || lmax < 0 || lmax ~= fix(lmax)
    error('csh_rot_wigner:lmax', 'lmax must be a non-negative integer scalar.');
end

if isnumeric(x) && isequal(size(x), [3 3])
    R = x;
    if norm(R.'*R - eye(3), 'fro') > 1e-10 || abs(det(R) - 1) > 1e-10
        error('csh_rot_wigner:R', ...
              'x must be a proper rotation matrix (orthogonal, det = +1).');
    end
    [alpha, beta, gamma] = mat2zyz(R);     % exact unitarity via unit-modulus phases
elseif isnumeric(x) && numel(x) == 3
    alpha = x(1); beta = x(2); gamma = x(3);
else
    error('csh_rot_wigner:x', ...
          'x must be a 3x3 rotation matrix or 3 ZYZ Euler angles.');
end

%% ---- band by band ----------------------------------------------------
M  = cell(1, lmax+1);
Ml = zeros(2*lmax+1, 2*lmax+1);
for l = 0:lmax
    D = wigD(l, alpha, beta, gamma);
    M{l+1} = D;
    idx = (lmax-l+1):(lmax+l+1);
    Ml(idx, idx) = D;
end

if nargout == 1
    M = Ml;                 % [Ml] = csh_rot_wigner(...) gives just the big matrix
end
end                                                  % ===== main =====


%% =====================================================================
function D = wigD(l, alpha, beta, gamma)
%WIGD  D^l_{m'm}(alpha,beta,gamma) = exp(-i m' alpha) d^l_{m'm}(beta) exp(-i m gamma)
m  = (-l:l).';
mp = (-l:l).';
D  = bsxfun(@times, bsxfun(@times, exp(-1i*mp*alpha), wigd(l, beta)), exp(-1i*m*gamma));
end


%% =====================================================================
function d = wigd(l, beta)
%WIGD  Reduced (small) Wigner d matrix d^l_{m'm}(beta), m', m = -l ... l.
%      Evaluated from the finite Jacobi-polynomial series; real, exact to
%      round-off.  Beta is reduced modulo 2*pi for accuracy.
beta = mod(beta, 2*pi);
if beta > pi, beta = beta - 2*pi; end      % now beta in (-pi, pi]
n  = 2*l + 1;
d  = zeros(n, n);
cb = cos(beta/2);
sb = sin(beta/2);

for i = 1:n
    mp = i - l - 1;                        % m' = -l ... l
    for j = 1:n
        m  = j - l - 1;                    % m  = -l ... l
        a1 = m - mp;                       % exponent of sin(beta/2)
        a2 = m + mp;                       % exponent of cos(beta/2)
        k  = l - m;                        % Jacobi degree

        % Normalisation (-1)^(m-m') sqrt((l+m)!(l-m)!/((l+m')!(l-m')!))
        N = (-1)^a1 * sqrt( (factrl(l+m)  * factrl(l-m)) / ...
                            (factrl(l+mp) * factrl(l-mp)) );

        % Jacobi P_k^{(a1,a2)}(x) at x = cos(beta), as a finite series.
        %   P_k^{(a,b)}(x) = sum_s (-1)^s C(k+a, k-s) C(k+b, s)
        %                    * ((1+x)/2)^(k-s) ((1-x)/2)^s
        % Written directly in terms of cos(beta/2), sin(beta/2) below, so
        % that the powers combine with the (cos(beta/2))^a2 (sin(beta/2))^a1
        % prefactor into non-negative integer exponents.
        tot = 0;
        for s = 0:k
            cjk = nchoosek(k + a1, k - s);     % zero when k+a1 < 0 handled below
            csk = nchoosek(k + a2, s);
            if k + a1 < 0 || k - s < 0 || k - s > k + a1, cjk = 0; end
            if k + a2 < 0 || s > k + a2,                  csk = 0; end
            cf  = cjk * csk;
            if cf == 0, continue, end
            e1  = 2*(k - s) + a2;              % exponent of cos(beta/2)
            e2  = 2*s       + a1;              % exponent of sin(beta/2)
            tot = tot + (-1)^s * cf * cb^e1 * sb^e2;
        end
        d(i, j) = N * tot;
    end
end

% small numerical cleanup: d is real and, at this sign convention, symmetric
d = real(d);
d = (d + d.') / 2;
end


%% =====================================================================
function f = factrl(n)
%FACTRL  Factorial of a non-negative integer, product form (avoids overflow
%        warnings for the O(30) values that occur here).
if n <= 1
    f = 1;
else
    f = prod(1:n);
end
end


%% =====================================================================
function [alpha, beta, gamma] = mat2zyz(R)
%MAT2ZYZ  ZYZ Euler angles from R = Rz(alpha) * Ry(beta) * Rz(gamma).
beta = atan2(sqrt(R(1,3)^2 + R(2,3)^2), R(3,3));
if sqrt(R(1,3)^2 + R(2,3)^2) > 1e-12
    alpha = atan2(R(2,3), R(1,3));
    gamma = atan2(R(3,2), -R(3,1));
else
    % gimbal lock (beta = 0 or pi): only alpha + gamma is determined
    alpha = 0;
    gamma = atan2(-R(1,2)*sign(R(3,3)), R(1,1)*sign(R(3,3)));
end
end
