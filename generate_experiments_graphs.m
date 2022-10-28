% This files generates each of the graphs used for the experiments shown in
%   the project report
n_nodes_values = [1024, 4096, 16384];
edges_per_node_values = [4, 8, 16];
seed_values = [1:1:5];

for n_nodes_index = 1:size(n_nodes_values, 2)
    n_nodes = n_nodes_values(n_nodes_index);
    fprintf("Nodes: %.0f\n", n_nodes);
    for edges_per_node_index = 1:size(edges_per_node_values, 2)
        edges_per_node = edges_per_node_values(edges_per_node_index);
        n_edges = n_nodes * edges_per_node;
        fprintf("\tEdges: %.0f\n", n_edges);
        for seed_index = 1:size(seed_values, 2)
            seed = seed_values(seed_index);
            fprintf("\t\tSeed: %.0f\n", seed);
            [A, b] = generate_problem_matrices(n_nodes, n_edges, seed, "skip", true);
        end
    end
end