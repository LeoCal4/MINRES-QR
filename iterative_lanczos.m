function [V, T, w_prime] = iterative_lanczos(A, V, T, w_prev, b, k, reorthogonalize, precon, D, C)
% A must be symmetric.
% Assuming V, T have been already created according to the algorithms
%   specifications. At step k = 1, both V and T are zero-matrices.
    %TODO add lucky breakdown check => x = 0
    if exist('reorthogonalize', 'var') == 0
       reorthogonalize = false;
    end
    if exist('precon', 'var') == 0
       precon = false;
    end
    if k == 1
        v_1 = b/norm(b); % O(n)
        V(:, 1) = v_1;
        w = A*v_1; % O(n^2)
        if precon ~= false
           w = multiply_preconditioner(w, D, C);
        end
        % alpha_1
        alpha_1 = v_1'*w;
        T(1, 1) = alpha_1; % O(n)
        w_prime = w - alpha_1*v_1; % O(n)
        % beta_1
        T(2, 1) = norm(w_prime); % O(n)
    else
        beta_k = T(k, k-1);
        T(k-1, k) = beta_k; 
        v_k = w_prev/beta_k; % O(n) 
        V(:, k) = v_k;
        % w_k
        w = A*v_k;
        if precon ~= false
           w = multiply_preconditioner(w, D, C);
        end
        w = w - beta_k*V(:, k-1); % O(n^2 + n)
        % alpha_k
        alpha_k = w'*v_k;
        T(k, k) = alpha_k; % O(n)
        w_prime = w - alpha_k*v_k; % O(n)
        if reorthogonalize
           w_prime = w_prime - V(:, 1:k) * (V(:, 1:k)' * w_prime);
        end
        T(k+1, k) = norm(w_prime); % O(n)
    end
end

