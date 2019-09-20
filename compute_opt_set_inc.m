function [ S_opt, omega, omega_list ] = compute_opt_set_inc( L_k, k, num_nodes_to_add, current_S_opt )
%   AUTHOR: Aamir Anis, USC
%   AUTHOR: Benjamin Girault, USC (assert, ordered list of sampled nodes, cleanup)
%   This function computes the optimal sampling set of a given size 
%   "S_opt_size" that maximizes the cutoff frequency.

% % %
% PARAMETER DESCRIPTION
% 
% INPUT
% L_k: kth power of Laplacian
% S_opt_size:  Desired size of the optimal set
% k: Power of Laplacian while computing cutoff, higher the order,
% greater the accuracy, but the complexity is also higher.
% 
% OUTPUT
% S_opt: Optimal set, ordered by selection
% omega: cutoff of the optimal set
% omega_list: List of computed cutoffs
% 
% % %

% Graph size
N = length(L_k);

% L_k Symmetric?
assert(ishermitian(L_k), 'StacUSC:ActiveSSLwithSampling:NonHermitian', 'Non Hermitian matrix L_k!');

% index vector
p = (1:N)';

% Initialization : If previous state available, initialize to that
if exist('current_S_opt','var')
    if numel(current_S_opt) == N
        current_S_opt = find(current_S_opt);
    end
    S_opt = current_S_opt;
else
    S_opt = [];
end
omega_list = zeros(num_nodes_to_add, 1);  

for iter = 1:num_nodes_to_add
    % create index vector for Sc from indicator functions
    S_complement = setdiff(p, S_opt);

    % compute minimum eigen-pair: efficient way
    [omega, y] = compute_cutoff(L_k, k, S_complement);

    % store a list of omega
    if nargout > 2
        omega_list(iter) = omega;
    end

    % find direction of maximum increment in reduced (|Sc|) dimensions
    [~,max_index] = max(abs(y));

    % Find corresponding node in N dimensions
    node_to_add = S_complement(max_index);

    % Update indicator function
    S_opt(end + 1) = node_to_add;  %#ok

    fprintf('Nodes added = %d...\n', numel(S_opt));
end

if nargout > 1
    omega = compute_cutoff(L_k, k, setdiff(p, S_opt));
    if nargout > 2
        omega_list = [omega_list; omega];
    end
end

end
