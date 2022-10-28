function [A, b] = generate_problem_matrices(...
        nodes, edges_parameter, seed, ...
        adj_matrix_filename, save_adj_matrix, ...
        verbose, min_max_D)
% Generate a random direct graph G and consequently the problem matrix A and
%   the b vector. 
% The graph is generated by creating the specified number of nodes and 
%   adding the specified number of edges (or arcs) between non-connected nodes,
%   making sure that each node has at least one connection.
%
% Inputs: 
%
% - nodes (int, scalar): the number of nodes, must be positive.
%
% - edges_parameters (real, scalar): defines the number of nodes in the graph. 
%   Can be setted to either one of these three options:
%   - edges_parameters > 1: number of edges;
%   - 0 < edges_parameters <= 1: density of the graph, defined as the
%       smallest nearest integer to density * (nodes^2 - nodes);
%   - -1 >= edges_parameters: the negative ratio of edges and nodes. For
%       example, setting nodes to 10 and edges_parameters to -5, the
%       resulting number of edges will be 10 * -(-5) = 50.
%   In any of these cases, the resulting number of edges must be greater
%   than the number of nodes.
%
% - seed (int, scalar, defaults to 0): the seed for RNG reproducibility.
% 
% Optional Inputs:
%
% - adj_matrix_filename (string, defaults to null): the filename (or relative path)
%       of the adjacency matrix to be opened. Specifying this avoids the
%       graph generation; instead, it creates G from the specified adjacency matrix
%       and extracts E to construct A.
%       If adj_matrix_filename is set to the empty string (""), the
%       filename is automatically set to
%       "matrices/adj_matrix_{n_nodes}_nodes_{n_edges}_edges_{seed}_seed.dat
%       If it is set to "skip", the value is discarded.
% 
% - save_adj_matrix (bool, defaults to false): if set to true, the
%       adjacency matrix of the generated graph is saved in a .dat file.
%       The naming convention follows that of the automatic name specified
%       for the adj_matrix_filename argument.
% 
% - verbose (bool, defaults to true): if set to true, the function prints
%       information regarding the matrices being created.
% 
% - min_max_D (real, [2 x 1 real vector], defaults to [1 2]): the range of
%       values to be used to draw random values from on the creation of D.
%
% Outputs:
%
% - A ([nodes+edges x nodes+edges] real sparse matrix): the main problem matrix, 
%       it is composed by [D E'; E 0], where D is a diagonal
%       positive definite matrix, while E is the incidence matrix of the graph.
%
% - b ([nodes+edges x 1] real column vector): defined as [c; d]. 
%       c ([edges x 1]) is simply created by generating #edges random values in range [0, 1].
%       To compley the constraint imposed by the linear dependency of E rows, 
%       the values of d ([nodes x 1]) must sum to 0; 
%       this is done by randomly generating the first half of d (d_1)
%       with positive values in the range [0, 1] and setting the second half 
%       to its opposite (d_2 = -d_1).
%

    % handle optional input variables
    if exist('seed', 'var') == 0
        seed = 0;
    end
    if exist('adj_matrix_filename', 'var') == 0
        adj_matrix_filename = "skip";
        open_adj_matrix = false;
    else
        open_adj_matrix = true;
    end
    if exist('save_adj_matrix', 'var') == 0
        save_adj_matrix = false;
    end
    if exist('verbose', 'var') == 0
        verbose = true;
    end
    if exist('min_max_D', 'var') == 0
        min_max_D = [1, 2];
    end
    
    % set seed
    rng(seed);
    
    % handle the edges_parameter value
    % if bigger than 1 it is number of edges
    if edges_parameter > 1
        edges = edges_parameter;
        density = edges / (nodes^2 - nodes);
    % if between 0 (not included) and 1 it is the density
    elseif edges_parameter > 0
        density = edges_parameter;
        edges = floor(density * (nodes^2 - nodes));
    % if smaller or equal to -1 it is -edge/nodes ratio, with edges > nodes
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
    
    % either read the graph from the file indicated by the provided filename
    %   or generate a new one
    if open_adj_matrix && strcmp(adj_matrix_filename, "skip") == 0
        [E, adj_matrix] = read_graph_from_file(adj_matrix_filename, nodes, edges, seed);
    else
        [E, adj_matrix] = generate_random_graph(nodes, edges);
    end
    
    % complete the creation of (sparse) A
    Z = sparse(nodes, nodes);
    A = [D E'; E Z];
    
    % randomly create b, enforcing the constraint of having its lower
    %   subvector d to sum to 0
    b = rand(nodes+edges, 1);
    half_d = b(edges+1:edges+nodes/2);
    b(edges+1+nodes/2:end) = -half_d;
    
    % save the adj_matrix
    if save_adj_matrix == true
        base_title = sprintf("adj_matrix_%i_nodes_%i_edges_%i_seed.dat", nodes, edges, seed);
        [i, j, val] = find(adj_matrix);
        writematrix([i j val; nodes nodes 0], base_title);
    end
end


function [E, adj_matrix, G] = read_graph_from_file(filename, nodes, edges, seed)
% Read the adjacency matrix from the relative path specified by the
%   filename variable. If filename is an empty string, the file is searched
%   in the matrices folder, using its default name specified by its nodes,
%   edges and seed.

    if filename == ""
       filename = sprintf("matrices/adj_matrix_%i_nodes_%i_edges_%i_seed.dat", nodes, edges, seed); 
    end
    fprintf("\tReading graph from %s\n", filename);
    adj_matrix_file = load(filename);
    adj_matrix = spconvert(adj_matrix_file);
    G = digraph(adj_matrix);
    E = sparse(G.incidence);
end


function [E, adj_matrix, G] = generate_random_graph(nodes, edges)
% Generates a random graph with #nodes and #edges specified by the input
%   arguments. The edges are randomly created by sampling two nodes at each
%   iteration and connecting them only if an edge between them is not
%   present. Each node has at least one outgoing edges, ensuring that the
%   resulting graph has only one connected component, making the rank of its
%   incidence matrix E equal to #nodes-1 (this is important for MINRES-QR convergence).

    % initialize E aka the incidence matrix
    E = sparse(nodes, edges);
    % adjacency matrix is created for the sole purpose of creating the
    %   graph object later
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
        % check if the edge already exist, without updating the col index
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
    %   edges for the first #nodes columns (this probably does not matter)
    random_E_cols_perm = randperm(edges);
    E = E(:, random_E_cols_perm);
    % generate the graph object
    G = digraph(adj_matrix);
end
