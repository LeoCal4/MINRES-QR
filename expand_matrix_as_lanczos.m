function [A_expanded] = expand_matrix_as_lanczos(A)
    % [A(1) A(2); A(2) rand(1, 1); 0 rand(1, 1)];
    A_size = size(A);
    A_expanded = [A zeros(A_size(1), 1); zeros(1, A_size(2)+1)];
    % alpha_{k}
    A_expanded(end-1, end) = rand(1, 1);
    % copy beta_{k} over alpha_{k}
    A_expanded(end-2, end) = A_expanded(end-1, end-1);
    % beta_{k+1}
    A_expanded(end, end) = rand(1, 1);
end

