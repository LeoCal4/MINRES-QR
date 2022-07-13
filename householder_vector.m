function [u, s] = householder_vector(x)
    s = norm(x);
    if x(1) >= 0
        s = -s;
    end
    v = x;
    v(1) = v(1) - s;
    u = v / norm(v);
end

