function [G] = generate_graph_matrix(size, density)
    if mod(size, 2) ~= 0
       disp("Size must be even")
       return
    end
    d = rand(size/2, 1);
    D = diag(d);
    E = sprand(size/2, size/2, density) > 0.5;
    Z = zeros(size/2, size/2);
    G = [D E'; E Z];
    G = full(G);
end
