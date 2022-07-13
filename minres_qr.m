function [x, residuals, setup_time, lanczos_times, qr_times, ...
    inversion_times, final_times, mat_mul1_times, mat_mul2_times] = minres_qr(A, b, optimize)
    % [culo, culo, setup_time, lanczos_times, qr_times, inv_times, final_times, , mat_mul1_times, mat_mul2_times]
    %tic
        size_A = size(A, 1);
        threshold = 1e-12;
        residuals = nan(1, size_A);
        lanczos_times = nan(1, size_A);
        qr_times = nan(1, size_A);
        inversion_times = nan(1, size_A);
        final_times = nan(1, size_A);
        % initialize Lanczos' matrices
        V = zeros(size_A, size_A+1);
        T = zeros(size_A+1, size_A);
        prev_w = 1;
        beta_1 = norm(b);
        % initialize QR matrices
        Q = eye(1);
        R = eye(1);
        setup_time = 0;
        mat_mul1_times = nan(size_A, 1);
        mat_mul2_times = nan(size_A, 1);
    %setup_time = toc;
    % compute algorithm iteration up to size of A
    for k = 1:size_A
        %disp("Iteration:")
        %k
        tic
            [V, T, prev_w] = iterative_lanczos(A, V, T, prev_w, b, k);
        lanczos_times(k) = toc;
        tic
            [Q, R, mat_mul1_times, mat_mul2_times] = iterative_QR(T(1:k+1, 1:k), Q, R, k, mat_mul1_times, mat_mul2_times, optimize);
        qr_times(k) = toc;
        %tic
            e_1 = zeros(k+1, 1);
            e_1(1) = 1;
            c = beta_1 * Q(:, 1:end-1)' * e_1;
            y_k = R(1:end-1, :) \ c;
        %inversion_times(k) = toc;
        %tic
            x = V(:, 1:k) * y_k;
            residual = norm(b - A*x);
            residuals(k) = residual;
        %final_times(k) = toc;
        if residual <= threshold
            break
        end
        %norm(T(1:k+1, 1:k) * y_k - beta_1 * e_1)
    end
    %norm(A*(b - A*x))
    fprintf("Completed in %g iterations\n", k)
end

