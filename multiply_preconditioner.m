function [b] = multiply_preconditioner(a, D, C, transpose_C)
    if exist('transpose_C', 'var') == 0
       transpose_C = false;
    end
    size_D = size(D, 1);
    b_1 = diag_inverter(D, a(1:size_D, :)); % D
    if transpose_C == true
       C = C';
    end
    b_2 = C \ a(size_D+1:end, :); % Sx = a => = -CC'x = a => x = -C'\C\a
    b = [b_1; b_2];
end
