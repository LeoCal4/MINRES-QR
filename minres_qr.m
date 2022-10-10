function [x, residuals, iter, residuals_precon] = minres_qr(A, b, reorthogonalize, precon, size_D)
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
        b = multiply_preconditioner(b, D_s, C, true);
        residuals_precon = nan(1, size(A, 1));
    end
    
    size_A = size(A, 1);
    threshold = 1e-9;
    residuals = nan(1, size_A);
    
    % initialize Lanczos' matrices
    V = zeros(size_A, size_A+1); % not sparse
    T = zeros(size_A+1, size_A); % sparse - prob not worth since most of the operations are just access and no mults
    prev_w = b;
    beta_1 = norm(b);
    
    % initialize QR matrices
    Q = eye(1); % not sparse
    R = eye(1); % sparse - needs to be assigned to create T_k in QR, so not worth to make it sparse
    c = zeros(1, 1);
    
    % iterate up to size of A at max
    for k = 1:size_A
%         fprintf("Iteration %i\n", k);
        
        % LANCZOS
        [V, T, prev_w] = iterative_lanczos(A, V, T, prev_w, k, reorthogonalize, precon, D_s, C);
        
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
    iter = k;
end
