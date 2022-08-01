function [G] = generate_graph_matrix(nodes, edges_parameter)
    % if bigger than 1 they are edges
    if edges_parameter > 1
        edges = edges_parameter;
        density = edges / (nodes^2 - nodes);
    % if between 0 (not included) and 1 it is the density
    elseif edges_parameter > 0
        density = edges_parameter;
        edges = floor(density * (nodes^2 - nodes));
    % if between smaller or equal to -1 it is -edge/nodes ratio, with
    %   edges > nodes
    elseif edges_parameter <= -1
        edges = nodes * -edges_parameter;
        density = edges / (nodes^2 - nodes);
    else
        error("Edges parameter must be != 0 and >= -1")
    end
    % check if there is at least one edge per node
    if edges < nodes 
        error("Not enough edges (nodes: %u - edges: %.0f)", nodes, edges)
    end
    fprintf("Generating graph with %u nodes and %.0f edges (density = %1.2f, edge/nodes = %1.2f)\n", nodes, edges, density, edges/nodes);
    fprintf("Final matrix size: %.0f x %.0f\n", nodes+edges, nodes+edges);
    d = rand(edges, 1);
    D = diag(d);
    E = zeros(nodes, edges);
    for col = 1:edges
        % ensure that each node has at least 1 (outgoing) edge
        if col < nodes +1
            first_random_row = col;
        % every node has at least one edge, so just take a random row
        else
            first_random_row = randi(nodes);
        end
        second_random_row = randi(nodes);
        % iterate until the -1 row is different from the 1 row
        while first_random_row == second_random_row
            second_random_row = randi(nodes);
        end
        E(first_random_row, col) = 1;
        E(second_random_row, col) = -1;
    end
    random_E_cols_perm = randperm(edges);
    E = E(:, random_E_cols_perm);
    Z = zeros(nodes, nodes);
    G = [D E'; E Z];
    G = full(G);
end
