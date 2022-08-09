function [Q_k, R_k] = iterative_QR(T_k, previous_Q, previous_R, k)
    % in any iteration > 1, the previous_Q is expanded for dimensional 
    %    consistency, and T_k is updated accordingly
    if k ~= 1
        % ===== OLD VERSION =====
        %previous_Q = [previous_Q zeros(k, 1)];
        %previous_Q = [previous_Q; zeros(1, k+1)];
        %previous_Q(end, end) = 1;
        %T_k = previous_Q' * T_k;
        % Optimization:
        %   can simplify since we already calculate the first k-1 rows and cols
        %   we can add the last row and cols to previous_R (R_k obtained from
        %   the previous iteration)
        % =======================
        previous_Q = [previous_Q; zeros(1, k)];
        new_col = previous_Q' * T_k(1:end, end); % O(n^2)
        Z = zeros(1, k-1);
        % new T_k
        T_k = [previous_R new_col; Z T_k(end, end)];
        previous_Q = [previous_Q zeros(k+1, 1)];
        previous_Q(end, end) = 1;
    end
    % calculate the current householder vector using only the last 2 elements
    %   of T_k k-th column
    [u_k, ] = householder_vector(T_k(k:end, k)); % O(1)
    % obtain the submatrix which will compose H_k
    H_k_prime = eye(2) - 2 * (u_k * u_k'); % O(1)
    % create H_k with padding blocks
    upper_right_block = zeros(k-1, 2);
    lower_left_block = zeros(2, k-1);
    H_k = [eye(k-1) upper_right_block; lower_left_block H_k_prime];
    % obtain R_k projecting T_k through H_k
    if k <= 2
        R_k = H_k * T_k; % O(1)
    else
        % we multiply just the lower right 2x2 matrix of T_k with
        % H_k_prime, since multiplying the whole matrices would only affect
        % said block
        T_k(end-1:end, end-1:end) = H_k_prime * T_k(end-1:end, end-1:end); % O(1)
        R_k = T_k;
    end
    % if it is the first iteration, H_k is already Q_k
    if k == 1
       Q_k = H_k;
       return
    end
    % otherwise, update the previous_Q with H_k 
    % ==== OLD VERSION ====
    %Q_k = previous_Q * H_k; 
    % Optimization:
    %   can simplfy by updating previous_Q and adding the last rows and
    %   cols. We multiply only the elements which are actually included in
    %   the whole matrix multiplication.
    % =====================
    Q_k = previous_Q;
    Q_k_col = Q_k(1:k, k); % saving the column to avoid overwriting it
    Q_k(1:k, k) = Q_k_col * H_k(k, k); % O(n)
    Q_k(1:k, k+1) = Q_k_col * H_k(k, k+1); % O(n)
    Q_k(k+1, k:k+1) = H_k(k+1, k:k+1);
end

