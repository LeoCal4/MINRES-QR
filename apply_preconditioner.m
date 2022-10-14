function [D_s, C, S, P] = apply_preconditioner(A, size_D, sparse, threshold)
    if exist('sparse', 'var') == 0
       sparse = false;
    end
    if exist('threshold', 'var') == 0
       threshold = 0;
    end
    D = A(1:size_D, 1:size_D);
    E = A(size_D+1:end, 1:size_D);
    fprintf("density(E): %2.5f\n", (nnz(E)*100)/(size(E, 1)*size(E, 2)));
    S = -E * (D \ E');
    fprintf("density(S): %2.5f\n", (nnz(S)*100)/numel(S));
    if sparse == true
        non_diag_mask = ~logical(eye(size(S)));
        S(non_diag_mask & S < threshold & S > -threshold) = 0;
        fprintf("post mask density(S): %2.5f\n", (nnz(S)*100)/numel(S));
    end
    D_s = sqrt(D);
    C = chol(-S);
    Z = zeros(size(C, 1), size_D);
    P = [D Z'; Z S];
end
