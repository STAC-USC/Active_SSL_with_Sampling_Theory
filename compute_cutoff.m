function [ omega, y ] = compute_cutoff( L_k, k, S_complement )
%COMPUTE_CUTOFF
%   AUTHOR: Aamir Anis, USC
%   AUTHOR: Benjamin Girault, USC (cleanup)
%   This function computes the cutoff frequency for a given
%   sampling set

% % %
% PARAMETER DESCRIPTION
% 
% INPUT
% L: kth power of Laplacian
% S_complement: Complement of the sampling set (list of non-sampled nodes indices)
% k: Power of Laplacian while computing cutoff, higher the order,
% greater the accuracy, but the complexity is also higher.
% 
% OUTPUT
% omega: cutoff frequency of the given set
% y: the associated eigenvector
% 
% % %

% L_k Symmetric?
assert(ishermitian(L_k), 'StacUSC:ActiveSSLwithSampling:NonHermitian', 'Non Hermitian matrix L_k!');

% compute minimum eigen-pair: efficient way
if nargout == 1
    omega = eigs(L_k(S_complement, S_complement), 1, 'sm');
else
    [y, omega] = eigs(L_k(S_complement, S_complement), 1, 'sm');
end
omega = abs(omega)^(1/k);

end
