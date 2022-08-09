function [x, residuals] = minres_qr(A, b, reorthogonalize, precon, size_D)
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
    % preconditing 
    D = 0;
    C = 0;
    if precon ~= false
       [D, C] = apply_preconditioner(A, size_D);
       b_1 = D \ b(1:size_D, :); % D
       b_2 = - C' \ (C \ b(size_D+1:end, :)); % S = -CC'
       b = [b_1; b_2];
    end
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
           y_1 = D \ y(1:size_D, 1);
           y_2 = - C' \ (C \ y(size_D+1:end, 1));
           y = [y_1; y_2];
        end
        residual = norm(b - y);
        residuals(k) = residual;
        if residual <= threshold
            break
        end
    end
    fprintf("Completed in %g iterations\n", k)
end

