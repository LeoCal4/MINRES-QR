function [A, b, E, D, G, A_t, b_t, E_t] = generate_graph_matrix(nodes, edges_parameter, seed, verbose, min_max_D)
    % handle optional variables
    if exist('seed', 'var') == 0
       seed = 0;
    end
    if exist('verbose', 'var') == 0
       verbose = true;
    end
    if exist('min_max_D', 'var') == 0
       min_max_D = [1, 2];
    end
    % set seed
    rng(seed);
    % if bigger than 1 they are edges
    if edges_parameter > 1
        edges = edges_parameter;
        density = edges / (nodes^2 - nodes);
    % if between 0 (not included) and 1 it is the density
    elseif edges_parameter > 0
        density = edges_parameter;
        edges = floor(density * (nodes^2 - nodes));
    % if smaller or equal to -1 it is -edge/nodes ratio, with
    %   edges > nodes
    elseif edges_parameter <= -1
        edges = nodes * -edges_parameter;
        density = edges / (nodes^2 - nodes);
    else
        error("Invalid edges parameter")
    end
    % check if there is at least one edge per node
    if edges < nodes - 1
        error("Not enough edges (nodes: %u - edges: %.0f)", nodes, edges)
    end
    if verbose
        fprintf("Generating graph with %u nodes and %.0f edges (density = %1.2f, edge/nodes = %1.2f)\n", nodes, edges, density, edges/nodes);
        fprintf("Final matrix size: %.0f x %.0f\n", nodes+edges, nodes+edges);
        fprintf("Values of D generated in range [%.0f, %.0f]\n",min_max_D(1), min_max_D(2));
    end
    % generate D according to its min and max values
    d = min_max_D(1) + (min_max_D(2)-min_max_D(1))*rand(edges, 1);
    D = sparse([1:1:edges], [1:1:edges], d);
    % initialize E aka the incidence matrix
    E = sparse(nodes, edges);
    % adjacency matrix is created for the sole purpose of creating the
    %   graph object later (which is not really needed)
    adj_matrix = zeros(nodes);
    % randomly create edges for the graph by iterating on the desired number 
    %   of edges, each time sampling two rows of the matrix 
    %   (hence two different nodes) and connecting them, i.e. setting one
    %   to 1 and the other to -1. There is no structure enforced, except
    %   for the fact that each node has at least one outgoing edge
    col = 1;
    while col < edges + 1
        % ensure that each node has at least 1 (outgoing) edge
        if col < nodes + 1
            first_random_row = col;
        % every node has at least one edge, so just take a random row
        else
            first_random_row = randi(nodes);
        end
        % iterate until the -1 row is different from the 1 row
        second_random_row = randi(nodes);
        while first_random_row == second_random_row
            second_random_row = randi(nodes);
        end
        % check if the edge already exist
        if adj_matrix(first_random_row, second_random_row) == 1
            continue
        end
        % update the matrices
        adj_matrix(first_random_row, second_random_row) = 1;
        E(first_random_row, col) = 1;
        E(second_random_row, col) = -1;
        col = col + 1;
    end
    % edge matrix columns are randomly permuted, to avoid having ordered
    %   edges for the first #nodes columns
    random_E_cols_perm = randperm(edges);
    E = E(:, random_E_cols_perm);
    % complete the creation of (sparse) A
    Z = sparse(nodes, nodes);
    A = [D E'; E Z];
    % generate the graph object
    G = digraph(adj_matrix);
    % randomly create b, enforcing the constraint of having its lower
    %   subvector to sum to 0
    b = rand(nodes+edges, 1);
    half_c = b(edges+1:edges+nodes/2);
    b(edges+1+nodes/2:end) = -half_c;
    % truncate the problem by one
    A_t = A(1:end-1, 1:end-1);
    b_t = b(1:end-1);
    E_t = E(1:end-1, :);
end
