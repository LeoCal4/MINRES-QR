function [D, C, S, P] = apply_preconditioner(A, size_D, sparse, threshold)
%function S = apply_preconditioner(A, size_D, sparse, threshold)
    if exist('sparse', 'var') == 0
       sparse = false;
    end
    if exist('threshold', 'var') == 0
       threshold = 0;
    end
    D = A(1:size_D, 1:size_D);
    E = A(size_D+1:end, 1:size_D);
    %fprintf("density(E): %2.5f\n", (nnz(E)*100)/(size(E, 1)*size(E, 2)));
    S = -E * (D \ E'); % a possible sparse approximation to do is to have nnz only where E has them
    %fprintf("density(S): %2.5f\n", (nnz(S)*100)/(size(S, 1)*size(S, 2)));
    if sparse == true
        E_mask = E(:, 1:size(E, 1)) ~= 0;
        S = S .* E_mask;
        %fprintf("post mask density(S): %2.5f\n", (nnz(S)/(size(S, 1)*size(S, 2)))*100);
        %S(S < threshold & S > -threshold) = 0;
    end
    C = chol(-S);
    Z = zeros(size(C, 1), size_D);
    P = [D Z'; Z S];
end
