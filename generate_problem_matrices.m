function [A_t, b_t, E] = generate_problem_matrices(...
        nodes, edges_parameter, ...
        seed, verbose, ...
        open_adj_matrix, adj_matrix_filename, save_adj_matrix, ...
        min_max_D)
    % handle optional variables
    if exist('seed', 'var') == 0
        seed = 0;
    end
    if exist('verbose', 'var') == 0
        verbose = true;
    end
    if exist("open_adj_matrix", "var") == 0
        open_adj_matrix = false;
    end
    if open_adj_matrix == true && exist('adj_matrix_filename', 'var') == 0
        % if this is empty, the adj matrix will be searched with its
        %   default name in the matrices folder, based on nodes, edges 
        %   and seed
        adj_matrix_filename = "";
    end
    if exist('save_adj_matrix', 'var') == 0
        save_adj_matrix = false;
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
        fprintf("Generating graph with %u nodes and %.0f edges (density = %1.2f, edge/nodes = %1.2f)\n", ...
            nodes, edges, density, edges/nodes);
        fprintf("Final matrix size: %.0f x %.0f\n", nodes+edges, nodes+edges);
        fprintf("Values of D generated in range [%.0f, %.0f]\n",min_max_D(1), min_max_D(2));
    end
    % generate D according to its min and max values
    d = min_max_D(1) + (min_max_D(2)-min_max_D(1))*rand(edges, 1);
    D = sparse([1:1:edges], [1:1:edges], d);
    if open_adj_matrix
        [E, adj_matrix, G] = read_graph_from_file(adj_matrix_filename, nodes, edges, seed);
    else
        [E, adj_matrix, G] = generate_random_graph(nodes, edges);
    end
    % complete the creation of (sparse) A
    Z = sparse(nodes, nodes);
    A = [D E'; E Z];
    % randomly create b, enforcing the constraint of having its lower
    %   subvector to sum to 0
    b = rand(nodes+edges, 1);
    half_c = b(edges+1:edges+nodes/2);
    b(edges+1+nodes/2:end) = -half_c;
    % truncate the problem by one
    A_t = A(1:end-1, 1:end-1);
    b_t = b(1:end-1);
    E_t = E(1:end-1, :);
    if save_adj_matrix == true
        base_title = sprintf("adj_matrix_%i_nodes_%i_edges_%i_seed.dat", nodes, edges, seed);
        [i, j, val] = find(adj_matrix);
        writematrix([i j val; nodes nodes 0], base_title);
    end
end

function [E, adj_matrix, G] = read_graph_from_file(filename, nodes, edges, seed)
    if filename == ""
       filename = sprintf("matrices/adj_matrix_%i_nodes_%i_edges_%i_seed.dat", nodes, edges, seed); 
    end
    fprintf("\tReading graph\n");
%     fileID = fopen(filename, "r");
%     formatSpec = '%i';
%     size_adj_matrix = [Inf Inf];
%     adj_matrix = fscanf(fileID,formatSpec, size_adj_matrix);
%     fclose(fileID);
    adj_matrix_file = load(filename);
    adj_matrix = spconvert(adj_matrix_file);
    G = digraph(adj_matrix);
    E = sparse(G.incidence);
end


function [E, adj_matrix, G] = generate_random_graph(nodes, edges)
    % initialize E aka the incidence matrix
    E = sparse(nodes, edges);
    % adjacency matrix is created for the sole purpose of creating the
    %   graph object later (which is not really needed)
    adj_matrix = sparse(nodes, nodes);
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
    % generate the graph object
    G = digraph(adj_matrix);
end
