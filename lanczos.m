function [V, T] = lanczos(A, b, k)
    
    if ~issymmetric(A)
        disp('A DEVE ESSERE SIMMETRICA OOOOOOOO')
        return
    end

    rows = size(A, 1);
    V = zeros(rows, k+1);
    T = zeros(k+1, k);
    %startup
    V(:, 1) = b/norm(b);
    T(1, 1) = V(:, 1)'*(A*V(:, 1));
    %T(1, 2) = norm(b);
    %variable for A*V(:, 1) to avoid double computation
    w_j_prime_prev = A*V(:, 1) - T(1, 1)*V(:, 1);
    T(2, 1) = norm(w_j_prime_prev);
    for j = 2 : k
        %do this here for 
        T(j-1, j) = T(j, j-1);
        %compute V(:, j): the jth vector of the orthonormal basis
        V(:, j) = w_j_prime_prev/T(j, j-1);
        w_j = A*V(:, j) - T(j, j-1)*V(:, j-1);
        alpha_j = w_j'*V(:, j);
        w_j_prime_prev  = w_j - alpha_j*V(:, j);
        %update T
        T(j+1, j) = norm(w_j_prime_prev);
        T(j, j) = alpha_j;
    end
    
    %T(end, end) = norm(w_j_prime_prev);
    %if T(end, end) ~= 0
    V(:, end) = w_j_prime_prev/T(end, end);
    %end
end

