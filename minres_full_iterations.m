n_nodes = 512;
edges_per_node_values = [8, 16];
% edges_per_node_values = [2, 3];
tol = 1e-09;
num_of_reps = 1;

for edges_per_node_index = 1:size(edges_per_node_values, 2)
    edges_per_node = edges_per_node_values(edges_per_node_index);
    n_edges = n_nodes * edges_per_node;
    maxit = n_nodes+n_edges;
    base_title = sprintf("minres_qr_%i_nodes_%i_edges", n_nodes, n_edges);
    fprintf(strcat(base_title, "\n\n"));
    iterations = nan(1, num_of_reps);
    times = nan(1, num_of_reps);
    residuals = nan(num_of_reps, maxit);
    for i = 1:num_of_reps
        [A, b, E, D, G, A_t, b_t, E_t] = generate_graph_matrix(n_nodes, n_edges, i);
        while rank(full(E)) ~= n_nodes - 1
            [A, b, E, D, G, A_t, b_t, E_t] = generate_graph_matrix(n_nodes, n_edges, i);
        end
        tic
        [x, res, iter] = minres_qr(A_t, b_t);
        times(1, i) = toc;
        iterations(1, i) = iter;
        residuals(i, 1:size(res, 2)) = res;
    end
    % write iterations
    writematrix(iterations, strcat(base_title, "_iterations.txt"));
    % write times
    writematrix(times, strcat(base_title, "_times.txt"));
    % write residuals
    writematrix(residuals, strcat(base_title, "_residuals.txt"));
    iter_means = sum(iterations, 2) / num_of_reps;
    fprintf("%i: %.3f\n", edges_per_node_index, iter_means);
end
