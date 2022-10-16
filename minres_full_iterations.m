n_nodes_values = [16384];
edges_per_node_values = [4];
% n_nodes_values = [1024, 4096, 16384];
% edges_per_node_values = [4, 8, 16];
tol = 1e-09;
num_of_reps = 5;

for n_nodes_index = 1:size(n_nodes_values, 2)
    n_nodes = n_nodes_values(n_nodes_index);
    fprintf("Nodes: %.0f\n", n_nodes);
    for edges_per_node_index = 1:size(edges_per_node_values, 2)
        edges_per_node = edges_per_node_values(edges_per_node_index);
        n_edges = n_nodes * edges_per_node;
        fprintf("\tEdges: %.0f\n", n_edges);
        maxit = n_nodes+n_edges;
        iterations = nan(1, num_of_reps);
        times = nan(1, num_of_reps);
        residuals = nan(num_of_reps, maxit);
        for i = 1:num_of_reps
            fprintf("\t\tRep: %.0f\n", i);
            [A_t, b_t] = generate_problem_matrices(n_nodes, n_edges, i, false, true);
            tic
            %%%%%% MINRES QR %%%%%% 
            [x, res, iter] = minres_qr(A_t, b_t, true, n_edges);
            %%%%%% MATLAB MINRES (PRECON %%%%%%%
%             [D_s, C] = apply_preconditioner(A_t, n_edges);
%             Z = sparse(n_nodes-1, n_edges);
%             M = [D_s Z'; Z C];
%             [x, flag, relres, iter, resvec] = minres(A_t, b_t, tol, 1000, M, M);
            times(1, i) = toc;
            iterations(1, i) = iter;
            residuals(i, 1:size(resvec, 1)) = resvec';
%             residuals(i, 1:size(res, 2)) = res;
        end
        base_title = sprintf("minres_qr_%i_nodes_%i_edges", n_nodes, n_edges);
        % write iterations
        writematrix(iterations, strcat(base_title, "_iterations.txt"));
        % write times
        writematrix(times, strcat(base_title, "_times.txt"));
        % write residuals
        writematrix(residuals, strcat(base_title, "_residuals.txt"));
%         iter_means = sum(iterations, 2) / num_of_reps;
%         fprintf("%i: %.3f\n", edges_per_node_index, iter_means);
    end
end