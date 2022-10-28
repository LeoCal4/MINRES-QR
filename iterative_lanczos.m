function [V, T, w_prime] = iterative_lanczos(A, V, T, w_prev, k, reorthogonalize, precon, D_s, C)
% Advances the Lanczos algorithm from iteration k-1 to iteration k,
%   provided the elements of the previous iteration.
%   At the first step k = 1, both V and T are zero-matrices and w_prev = b.
%   The optional arguments regulate both the possibilty to apply
%   reorthogonalization and preconditioning to the process.
%
% Inputs:
%
% - A ([(nodes-1)+edges x (nodes-1)+edges] real sparse matrix): the main problem matrix.
%
% - V ([(nodes-1)+edges x k-1] real matrix or scalar): the orthogonal matrix V_{k-1}
%       containing the vectors spanning the Krylov subspace up to dimension k-1.
%       During the first iteration, this is just 0.
%
% - T ([k x k-1] real matrix or scalar): the tridiagonal matrix T_{k-1}
%       During the first iteration, this is just 0.
%
% - w_prev ([(nodes-1)+edges x 1] real vector): the vector w_{k-1} obtained during the previous
%       iteration. It is needed to create v_k. 
%       IMPORTANT: During the first iteration, this is the vector b.
%
% - k (int, scalar): the current iteration number.
%
% Optional inputs:
%
% - reorthogonalize (bool, defaults to false): whether to reorthogonalize
%       the system or not. This noticebly increases the cost of the steps,
%       while making the algorithm backward stable.
%
% - precon (bool, defaults to false): 
%
% - D_s ([edges x edges] real matrix, defaults to 0): the square root of D.
%
% - C ([nodes-1 x nodes-1] real matrix, defaults to 0): the Cholesky factor of -S.
%
% Outputs:
%
% - V ([(nodes-1)+edges x k] real matrix or scalar): the updated matrix V_k.
%
% - T ([k+1 x k] real matrix or scalar): the updated tridiagonal matrix T_k.
%
% - w_prime ([(nodes-1)+edges x 1] real vector): the current w'; it becomes the next w_prev.
%
    
    % optional variables check
    if exist('reorthogonalize', 'var') == 0
       reorthogonalize = false;
    end
    if exist('precon', 'var') == 0
       precon = false;
    end
    if exist('D_s', 'var') == 0
       D_s = 0;
    end
    if exist('C', 'var') == 0
       C = 0;
    end
    
    % calculate v_k and update T_k and V_k matrices
    if k ~= 1
        beta_k = T(k, k-1);
        T(k-1, k) = beta_k;
        v_k = w_prev/beta_k;
    else
        % w_prev = b at k = 1, hence v_1 = b / norm(b)
        v_k = w_prev/norm(w_prev);
    end
    V(:, k) = v_k;
    
    % calculate w_k
    if precon
        w = multiply_preconditioner(v_k, D_s, C);      % B^-1 * v_k
        w = A * w;  % sparse multiplication            % A * B^-1 * v_k
        w = multiply_preconditioner(w, D_s, C, true);  % B^-T * A * B^-1 * v_k
    else
        w = A * v_k; % sparse multiplication
    end

    if k ~= 1
        w = w - beta_k*V(:, k-1);
    end
    
    % calculate \alpha_k and w' (which will become the next w_prev)
    alpha_k = w'*v_k;
    T(k, k) = alpha_k; 
    w_prime = w - alpha_k*v_k;
    
    % apply reorthogonalization if specified
    if reorthogonalize && k > 1
       w_prime = w_prime - V(:, 1:k) * (V(:, 1:k)' * w_prime); % O(n - 2kn+ 2kn)
    end
    
    T(k+1, k) = norm(w_prime);
    
end

