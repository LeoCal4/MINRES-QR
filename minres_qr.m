function [x, residuals, residuals_precon] = minres_qr(A, b, reorthogonalize, precon, size_D)
    if exist('reorthogonalize', 'var') == 0
       reorthogonalize = false;
    end
    if exist('precon', 'var') == 0
       precon = false;
    end
    if ~is_symm(A)
        disp('A DEVE ESSERE SIMMETRICA OOOOOOOO')
        return
    end
    
    % saving the original A and b
    b_start = b;
    A_start = A;
    
    % preconditioning 
    D_s = 0;
    C = 0;
    if precon == true
        A = A(1:end-1, 1:end-1);
        b = b(1:end-1, :);
        [D_s, C] = apply_preconditioner(A, size_D);
        %b_start = b;
        b = multiply_preconditioner(b, D_s, C, true);
        residuals_precon = nan(1, size(A, 1));
        
        % TEST PRECON
%         B = zeros(size(D_s, 1) + size(C, 1));
%         B(1:size(D_s, 1), 1:size(D_s, 1)) = D_s;
%         B(size(D_s, 1)+1:end, size(D_s, 1)+1:end) = C;
%         A_precon = B' \ A / B;
    end
    
    size_A = size(A, 1);
    threshold = 1e-9;
    residuals = nan(1, size_A);
    
    % initialize Lanczos' matrices
    V = zeros(size_A, size_A+1);
    T = zeros(size_A+1, size_A);
    prev_w = b;
    beta_1 = norm(b);
    
    % initialize QR matrices
    Q = eye(1);
    R = eye(1);
    c = zeros(1, 1);
    
    % iterate up to size of A at max
    for k = 1:size_A
%         fprintf("Iteration %i\n", k);
        
        % LANCZOS
        [V, T, prev_w] = iterative_lanczos(A, V, T, prev_w, k, reorthogonalize, precon, D_s, C);
%         [V, T, prev_w] = iterative_lanczos_inner(A, V, T, prev_w, k, D_s, C, false);
        
%         % TEST LANCZOS
%         [V_1, T_1, prev_w_1] = iterative_lanczos(A, V, T, prev_w, k, reorthogonalize, precon, D_s, C);
%         [V_2, T_2, prev_w_2] = iterative_lanczos(A_precon, V, T, prev_w, k);
%         fprintf("Difference between mult_prec and direct Lanczos: %f\n", norm(prev_w_1 - prev_w_2)/norm(prev_w_2));
%         V = V_2;
%         T = T_2;
%         prev_w = prev_w_2;
        
        % QR
        [Q, R] = iterative_QR(T(1:k+1, 1:k), Q, R, k);
        
        % BACKSUBSTITUTION
        %  use the reduced version of the QR, using just Q_1 and R_1
        %      from Q = [Q_1 Q_2] and R = [R_1; 0]
        c(k, 1) = beta_1 * Q(1, k);
        y_k = R(1:end-1, :) \ c; % O(n)
        
        % SOLUTION
        x_k = V(:, 1:k) * y_k;
        
        if precon == true
            % solution of the original truncated system
            x_k = multiply_preconditioner(x_k, D_s, C); % x = B^-1 * \hat(x)
            y_p = A * x_k; % y_p = A * B^-1 * \hat(x)
            % A*x of the preconditioned system
            y_p = multiply_preconditioner(y_p, D_s, C, true); % y_p = B^-T * A * B^-1 * \hat(x)
            residuals_precon(k) = norm(b - y_p);
%             x = multiply_preconditioner(x_k, D_s, C);
            x_k(end+1) = 0; % solution of the original system
        end
        
        % A*x of the original system
        y = A_start * x_k;
        residual = norm(b_start - y);
        residuals(k) = residual;
        if precon == false && residual <= threshold
            break
        elseif precon == true && residuals_precon(k) <= threshold
            residuals_precon = residuals_precon(1:k);
            break
        end
    end
    x = x_k;
    residuals = residuals(1:k);
    fprintf("Completed in %g iterations\n", k)
end


function [V, T, w_prime] = iterative_lanczos_inner(A, V, T, w_prev, k, D_s, C, full_precon)
    % calculate v_k
    if k ~= 1
        beta_k = T(k, k-1);
        T(k-1, k) = beta_k;
        v_k = w_prev/beta_k; % O(n)
    else
         v_k = w_prev/norm(w_prev); % w_prev = b at k = 1
    end
    V(:, k) = v_k;
    
    % calculate w_k
    if full_precon == true % A is already B^-T * A * B^-1
        w = A * v_k;
    else
        w = multiply_preconditioner(v_k, D_s, C);    % B^-1 * v_k
        w = A * w;                                   % A * B^-1 * v_k
        w = multiply_preconditioner(w, D_s, C, true);  % B^-T * A * B^-1 * v_k
    end
    % calculate w_k
%     w_direct = A_precon * v_k;
%     w = multiply_preconditioner(v_k, D_s, C);    % B^-1 * v_k
%     w = A * w;                                   % A * B^-1 * v_k
%     w = multiply_preconditioner(w, D_s, C, true);  % B^-T * A * B^-1 * v_k
%     fprintf("w direct - w mp: %.10f\n", norm(w_direct - w)/norm(w_direct));

    if k ~= 1
        w = w - beta_k*V(:, k-1); % O(n^2 + n)
    end
    
    alpha_k = w'*v_k; % O(n^2)
    T(k, k) = alpha_k; 
    w_prime = w - alpha_k*v_k; % O(2n)
    T(k+1, k) = norm(w_prime); % O(n)
end
