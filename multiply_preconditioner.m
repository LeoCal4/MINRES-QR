function [b] = multiply_preconditioner(a, D, C)
    size_D = size(D, 1);
    b_1 = D \ a(1:size_D, :); % D
    b_2 = - (C' \ (C \ a(size_D+1:end, :))); % Sx = a => = -CC'x = a => x = -C'\C\a
    b = [b_1; b_2];
    %S = -C'*C;
    %Z = zeros(size(C, 1), size_D);
    %M = [D Z'; Z S];
    %a = M \ a;
end
