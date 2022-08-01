function [V, T, w_prime] = iterative_lanczos(A, V, T, w_prev, b, k, reorthogonalize)
%A must be symmetric.
%Assuming V, T have been already created according to the algorithms
%specifications. At step k = 1, both V and T are zero-matrices.
    %TODO add lucky breakdown check => x = 0
    if exist('reorthogonalize', 'var') == 0
       reorthogonalize = false;
    end
    if k == 1
        % How to pass b? We only need it at the first iteration
        v_1 = b/norm(b); % O(n)
        V(:, 1) = v_1;
        w = A*v_1; % O(n^2)
        % alpha_1
        alpha_1 = v_1'*w;
        T(1, 1) = alpha_1; % O(n)
        w_prime = w - alpha_1*v_1; % O(n)
        % beta_1
        T(2, 1) = norm(w_prime); % O(n)
    else
      %   p1      = Operator * v1  -  beta1 * v0 => w = A * v_1 - beta_1 * v_0
      %   alpha1  = v1'p1 => T(1, 1) = a_1 = v_1 * w'
      %   q2      = p2  -  alpha1 * v1 => w_prime = w - T(1, 1)[aka a_1] * v_1
      %   beta2^2 = q2'q2 -
      %   v2      = (1/beta2) q2 => norm(w_prime)
        beta_k = T(k, k-1);
        T(k-1, k) = beta_k; 
        v_k = w_prev/beta_k; % O(n) 
        V(:, k) = v_k;
        % w_k
        w = A*v_k - beta_k*V(:, k-1); % O(n^2 + n)
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

