function [V, T, w_prime] = iterative_lanczos(A, V, T, w_prev, k, reorthogonalize, precon, D_s, C)
% A must be symmetric.
% Assuming V, T have been already created according to the algorithms
%   specifications. At step k = 1, both V and T are zero-matrices and w_prev = b.

    %TODO add lucky breakdown check => x = 0
    
    % optional variables check
    if exist('reorthogonalize', 'var') == 0
       reorthogonalize = false;
    end
    if exist('precon', 'var') == 0
       precon = false;
    end
    
    % calculate v_k
    if k ~= 1
        beta_k = T(k, k-1);
        T(k-1, k) = beta_k;
        v_k = w_prev/beta_k; % O(n)
    else
        % w_prev = b at k = 1
        v_k = w_prev/norm(w_prev); % O(n)
    end
    V(:, k) = v_k;
    
    % calculate w_k
    if precon == true % O(2m + t^2 + n^2)
        w = multiply_preconditioner(v_k, D_s, C);      % B^-1 * v_k
        w = A * w;                                     % A * B^-1 * v_k
        w = multiply_preconditioner(w, D_s, C, true);  % B^-T * A * B^-1 * v_k
    else
        w = A * v_k; % O(n^2)
    end

    if k ~= 1
        w = w - beta_k*V(:, k-1); % O(2n)
    end
    
    alpha_k = w'*v_k; % O(2n)
    T(k, k) = alpha_k; 
    w_prime = w - alpha_k*v_k; % O(2n)
    if reorthogonalize && k > 1
       w_prime = w_prime - V(:, 1:k) * (V(:, 1:k)' * w_prime); % O(n - 2kn+ 2kn)
    end
    T(k+1, k) = norm(w_prime); % O(n)

% function [V, T, w_prime] = iterative_lanczos(A, V, T, w_prev, b, k, reorthogonalize, precon, D, C)
%     if k == 1
%         v_1 = b/norm(b); % O(n)
%         V(:, 1) = v_1;
%         if precon ~= false
%            w = multiply_preconditioner(v_1, D, C);
%            %v_1(size(D, 1)+1:end) = -v_1(size(D, 1)+1:end);
%         end
%         w = A*w; % O(n^2)
%         if precon ~= false
%            w = multiply_preconditioner(w, D, C, true);
%         end
%         % alpha_1
%         alpha_1 = v_1'*w;
%         T(1, 1) = alpha_1; % O(n)
%         w_prime = w - alpha_1*v_1; % O(n)
%         % beta_1
%         T(2, 1) = norm(w_prime); % O(n)
%     else
%         beta_k = T(k, k-1);
%         T(k-1, k) = beta_k; 
%         v_k = w_prev/beta_k; % O(n) 
%         V(:, k) = v_k;
%         % w_k
%         if precon ~= false
%            w = multiply_preconditioner(v_k, D, C);
%            %v_k(size(D, 1)+1:end) = -v_k(size(D, 1)+1:end);
%         end
%         w = A*w; % O(n^2)
%         if precon ~= false
%            w = multiply_preconditioner(w, D, C, true);
%         end
%         w = w - beta_k*V(:, k-1); % O(n^2 + n)
%         alpha_k = w'*v_k; % O(n^2)
%         T(k, k) = alpha_k; 
%         w_prime = w - alpha_k*v_k; % O(2n)
%         if reorthogonalize
%            w_prime = w_prime - V(:, 1:k) * (V(:, 1:k)' * w_prime);
%         end
%         T(k+1, k) = norm(w_prime); % O(n)
%     end
end

