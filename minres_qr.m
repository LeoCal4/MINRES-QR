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
    b_start = b;
    % preconditioning 
    D = 0;
    C = 0;
    if precon == true
       A = A(1:end-1, 1:end-1);
       b = b(1:end-1, :);
       [D, C] = apply_preconditioner(A, size_D);
       %b = multiply_preconditioner(b, D, C);
       residuals_precon = nan(1, size(A, 1));
    end
    size_A = size(A, 1);
    threshold = 1e-9;
    residuals = nan(1, size_A);
    % initialize Lanczos' matrices
    V = zeros(size_A, size_A+1);
    T = zeros(size_A+1, size_A);
    prev_w = 1;
    beta_1 = norm(b);
    % initialize QR matrices
    Q = eye(1);
    R = eye(1);
    c = zeros(1, 1);
    % compute algorithm iteration up to size of A
    for k = 1:size_A
        [V, T, prev_w] = iterative_lanczos(A, V, T, prev_w, b, k, reorthogonalize, precon, D, C);
        [Q, R] = iterative_QR(T(1:k+1, 1:k), Q, R, k);
        % use the reduced version of the QR, using just Q_1 and R_1
        % from Q = [Q_1 Q_2] and R = [R_1; 0]
        c(k, 1) = beta_1 * Q(1, k);
        y_k = R(1:end-1, :) \ c; % O(n)
        x = V(:, 1:k) * y_k;
        y = A*x;
        if precon ~= false
           y_p = multiply_preconditioner(y, D, C);
           residuals_precon(k) = norm(b - y_p);
           y(end+1) = 0;
           x(end+1) = 0;
        end
        residual = norm(b_start - y);
        residuals(k) = residual;
        if residual <= threshold
            residuals = residuals(1:k);
            break
        end
    end
    fprintf("Completed in %g iterations\n", k)
end

