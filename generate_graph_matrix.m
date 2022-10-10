function [A, b, E, D, G, A_t, b_t, E_t] = generate_graph_matrix(nodes, edges_parameter, seed, min_max_D, polytree)
    if exist('seed', 'var') == 0
       seed = 0;
    end
    if exist('min_max_D', 'var') == 0
       min_max_D = [1, 2];
    end
    if exist('polytree', 'var') == 0
       polytree = false;
    end
    rng(seed);
    %if polytree == true
    %   edges = nodes - 1; 
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
    if edges < nodes - 1
        error("Not enough edges (nodes: %u - edges: %.0f)", nodes, edges)
    end
    fprintf("Generating graph with %u nodes and %.0f edges (density = %1.2f, edge/nodes = %1.2f)\n", nodes, edges, density, edges/nodes);
    fprintf("Final matrix size: %.0f x %.0f\n", nodes+edges, nodes+edges);
    fprintf("Values of D generated in range [%.0f, %.0f]\n",min_max_D(1), min_max_D(2));    
    d = min_max_D(1) + (min_max_D(2)-min_max_D(1))*rand(edges, 1);
%     D = sparse(diag(d));
    D = sparse([1:1:edges], [1:1:edges], d);
    E = sparse(nodes, edges);
    adj_matrix = zeros(nodes);
    % TODO add check to not add an edge which is already present
    for col = 1:edges
        % ensure that each node has at least 1 (outgoing) edge
        if col < nodes + 1
            first_random_row = col;
        % every node has at least one edge, so just take a random row
        else
            first_random_row = randi(nodes);
        end
        
        % avoid having outgoing edges from the last node in case of polytree
        if polytree == true && first_random_row == nodes
            while first_random_row == nodes
                first_random_row = randi(nodes);  
            end
        end
        
        if polytree == false
            second_random_row = randi(nodes);
        else
            % restrict selection to nodes which come after the first one
            second_random_row = randi([first_random_row nodes]);
        end
        % iterate until the -1 row is different from the 1 row
        while first_random_row == second_random_row
            if polytree == false
                second_random_row = randi(nodes);
            else
                second_random_row = randi([first_random_row nodes]);
            end
            % fprintf("Node1: %i - Node2: %i\n", first_random_row, second_random_row);
        end
        adj_matrix(first_random_row, second_random_row) = 1;
        E(first_random_row, col) = 1;
        E(second_random_row, col) = -1;
    end
    random_E_cols_perm = randperm(edges);
    E = E(:, random_E_cols_perm);
    Z = zeros(nodes, nodes);
    A = [D E'; E Z];
%     A = full(A);
    G = digraph(adj_matrix);
    b = rand(nodes+edges, 1);
    half_c = b(edges+1:edges+nodes/2);
    b(edges+1+nodes/2:end) = -half_c;
    A_t = A(1:end-1, 1:end-1);
    b_t = b(1:end-1);
    E_t = E(1:end-1, :);
end
