function [C] = matrix_equality(A,B)
    C = abs(A-B) < 1e4 * eps(min(abs(A), abs(B)))
end

