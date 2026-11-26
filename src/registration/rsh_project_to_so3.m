function R_proj = rsh_project_to_so3(M)
%RSH_PROJECT_TO_SO3  Nearest special orthogonal matrix (Frobenius norm).
%
%   R_proj = RSH_PROJECT_TO_SO3(M)
%
% Projects a (2*l+1) x (2*l+1) matrix M to the nearest orthogonal matrix
% with determinant +1 in Frobenius norm.
% This is generally not guaranteed to yield a *valid RSH rotation matrix*
% (i.e., one that is part of the specific irreducible representation).
% However, for standard numerical stability, it's used to ensure
% orthogonality and unit determinant.
%
% INPUTS:
%   M : (2*l+1) x (2*l+1) matrix
%
% OUTPUTS:
%   R_proj : Nearest (2*l+1) x (2*l+1) orthogonal matrix with det=+1

    [U, ~, V] = svd(M);
    
    % Ensure determinant is +1.
    % The projection to nearest orthogonal matrix is U*V'.
    % If det(U*V') is -1, then flip the sign of one column of U (or V) to make it +1.
    if det(U * V.') < 0
        U(:, end) = -U(:, end);
    end
    
    R_proj = U * V.';
    
    % Note: While this projects to SO(N), it does NOT guarantee that the resulting
    % matrix is a valid SO(3) rotation matrix *in the RSH basis*.
    % For that, it should be of the form exp(W_ell).
    % However, if the initial R0 is a valid rotation and all steps are exp(W_ell),
    % this is mainly a numerical safeguard and *shouldn't* be hit often or significantly.
end
