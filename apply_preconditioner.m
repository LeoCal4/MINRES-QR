function [D_s, C, S] = apply_preconditioner(A, size_D, sparse_prec, threshold, percentage_threshold)
    if exist('sparse_prec', 'var') == 0
       sparse_prec = false;
    end
    if exist('threshold', 'var') == 0
       threshold = 0;
    end
    if exist('percentage_threshold', 'var') == 0
       percentage_threshold = 0;
    end
    D = A(1:size_D, 1:size_D);
    E = A(size_D+1:end, 1:size_D);
    fprintf("density(E): %.5f\n", (nnz(E)*100)/numel(E));
    S = -E * (D \ E'); %fast
    fprintf("density(S): %.5f\n", (nnz(S)*100)/numel(S));
    if sparse_prec == true
        non_diag_mask = ~logical(eye(size(S)));
        % S seems to have negative values only on diagonal, hence the
        %   threshold is just for positive values
%         fprintf("Positive S mean: %.5f - negative S mean: %.5f\n", ...
%                 mean(mean(full(S(S > 0 & non_diag_mask)))), mean(mean(full(S(S < 0 & non_diag_mask)))));
%         fprintf("Positive S min-max: [%.5f, %.5f]\n", ...
%                 min(min(full(S(S > 0 & non_diag_mask)))), max(max(full(S(S > 0 & non_diag_mask)))));
        full_masked_S = full(S(S > 0 & non_diag_mask));
        min_value = min(min(full_masked_S));
        threshold = min_value + ((max(max(full_masked_S)) - min_value)/100)*percentage_threshold; %fast
        fprintf("Threshold: %d\n", threshold);
        S(non_diag_mask & S < threshold) = 0; % slow
        fprintf("post mask density(S): %2.5f\n", (nnz(S)*100)/numel(S));
    end
    D_s = sqrt(D);
    C = chol(-S); % slow
    fprintf("density(C): %2.5f\n", (nnz(C)*100)/numel(C));
%     Z = sparse(size(C, 1), size_D);
%     P = [D Z'; Z S];
end
