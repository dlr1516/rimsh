function [V, D] = eig_sort_imag(M)
% EIG_SORT_IMAG
% Computes the eigenvalues and eigenvectors of an input matrix and returns 
% them sorted accordig to their imaginary value. 
% It has been designed for Spherical Harmonics Rotation matrices whose
% eigenvalues (in the case of matrices of order l) are:
% -l * i, (-l+1) * i, ..., -2 * i, -i, 0, i, 2* i, ..., l * i   
%
% INPUTS:
% - M: the matrix whose eigenvalues are to be computed
%
% OUTPUTS:
% - V: a full matrix V whose columns are the corresponding eigenvectors 
%      sorted by the increasing value of imaginary part of eigenvalues
% - D: diagonal matrix D of eigenvalues sorted by the increasing values of 
%      their imaginary part
%
    [V, D] = eig(M);
    [~,idx] = sort(imag(diag(D)));
    P = eye(size(D));
    P = P(idx, :);
    D = P * D * P';
    V = V * P';
end