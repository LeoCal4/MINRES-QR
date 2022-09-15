function [A, E, D, G] = read_graph_from_file()
    fprintf("Reading graph\n");
    fileID = fopen("adj_matrix_50_005.txt", "r");
    formatSpec = '%i';
    size_adj_matrix = [50 Inf];
    adj_matrix = fscanf(fileID,formatSpec, size_adj_matrix);
    fclose(fileID);
    G = digraph(adj_matrix);
    E = G.incidence;
    n_nodes = size(E, 1);
    n_edges = size(E, 2);
    d = rand(n_edges, 1);
    D = diag(d);
    Z = zeros(n_nodes);
    A = [D E'; E Z];
end

