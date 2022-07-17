function [prec_A, prec_b] = apply_preconditioner(A, b, optimize, sparse, threshold)
    if exist('sparse', 'var') == 0
       sparse = false;
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
    %eig(inv(P)*A)
    [L, ] = ldl(P);
    % precomputing the inverse is faster than reinverting at each operation
    if optimize == false
        prec_A = L \ A / L';
        prec_b = L \ b;
    else
        inv_L = inv(L);
        prec_A = inv_L * A * inv_L';
        prec_b = inv_L * b;
    end
end
