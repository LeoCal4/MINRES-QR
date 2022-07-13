function [is_symmetric] = is_symm(A)
    for i = 1:size(A, 1)
        for j = i+1:size(A, 1)
            if abs(A(i, j) - A(j, i)) > 1e-13
                abs(A(i, j) - A(j, i))
                is_symmetric = false;
                return
            end
    end
    is_symmetric = true;
end

