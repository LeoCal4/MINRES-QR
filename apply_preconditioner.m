function [prec_A, prec_b] = apply_preconditioner(A, b, optimize, sparse, threshold)
    if exist('optimize', 'var') == 0
       optimize = true;
    end
    if exist('sparse', 'var') == 0
       sparse = false;
    end
    if exist('threshold', 'var') == 0
       threshold = 0;
    end
    size_A = size(A, 1);
    if mod(size_A, 2) ~= 0
        disp("Matrix size is not even")
        return 
    end
    D = A(1:size_A/2, 1:size_A/2);
    E = A((size_A/2)+1:end, 1:size_A/2);
    Z = zeros(size_A/2);
    if optimize == false
        S = -E * inv(D) * E';
    else
        S = -E * (D \ E');
    end
    if sparse == true
        S(S < threshold & S > -threshold) = 0;
    end
    P = [D Z; Z S];
    %eig(P\A)
    % using decomposition() function is faster than precomputing the inverse 
    % as the decomposition does not explicitely build L and D (probably)
    [L, ] = ldl(P);
    if optimize == false
        inv_L = inv(L);
        prec_A = inv_L * A * inv_L';
        prec_b = inv_L * b;
    else
        dL = decomposition(L);
        prec_A = dL \ A / dL';
        prec_b = dL \ b;
    end
end
