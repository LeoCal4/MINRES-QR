function [D_s, C, P] = create_preconditioner(A, n_edges, threshold, build_full_precon, verbose)
% Creates the Schur Complement preconditioner and its components, 
%   needed to multiply it efficiently and indirectly, while maintaining symmetry. 
%   The preconditioner is defined as [D 0; 0 S], where D is the diagonal
%   positive matrix from the problem matrix A, while S is its Schur
%   complement with respect to D, defined as S = - E D^-1 E. If E (the
%   incidence matrix) has full rank (#nodes - 1), the matrix S is negative
%   definite. This needs to be true, as otherwise Cholesky cannot be
%   applied to -S.
%   The two blocks composing the matrix are factorized using Cholesky, 
%   which is applied to D (indirectly) to obtain D_s, i.e. square root of D,
%   and to -S (the minus sign is caused by the negative positiveness) to obtain C.
%
% Inputs:
%
% - A ([nodes+edges x nodes+edges] real matrix): the main problem matrix, 
%       composed by [D E'; E 0]
% 
% - n_edges (int, scalar): the number of edges in the problem graph. This
%       is needed to know how to extract D from A.
%
% Optional inputs:
%
% - threshold (real, scalar, defaults to -1): specifies the desired non-negative threshold for non-diagonal 
%       elements in the S block.
%       It can be defined in the two following ways:
%       - threshold > 0: the desired threshold value;
%       - threshold = 0: the threshold value is set to be the mean of the
%           non-diagonal values of S.
%
% - build_full_precon (bool, defaults to false): whether to build the full
%       preconditioner matrix P.
%
% - verbose (bool, defaults to false): prints information about the
%       matrices densities and S values range.
%
% Outputs:
%
% - D_s ([edges x edges] real matrix): the square root of D.
% 
% - C ([nodes x nodes] real matrix): the Cholesky factor of -S.
% 
% - P ([nodes+edges x nodes+edges], real matrix, optional output): 
%       the full preconditioner matrix.
%

    % handle optional input arguments
    if exist('threshold', 'var') == 0
        threshold = -1;
    end
    if exist('build_full_precon', 'var') == 0
       build_full_precon = false; 
    end
    if exist('verbose', 'var') == 0
       verbose = false; 
    end
    D = A(1:n_edges, 1:n_edges);
    E = A(n_edges+1:end, 1:n_edges);
    S = -E * (D \ E');
    if verbose
        fprintf("density(E): %.5f\n", (nnz(E)*100)/numel(E));
        fprintf("density(S): %.5f\n", (nnz(S)*100)/numel(S));
    end
    
    % sparsify S using threshold
    if threshold >= 0
        % the diagonal elements of S cannot be changed, as otherwise the
        %   matrix would not be neg def anymore
        non_diag_mask = ~logical(eye(size(S)));
        if verbose
            fprintf("Positive S mean: %.5f - negative S mean: %.5f\n", ...
                    mean(mean(full(S(S > 0 & non_diag_mask)))), mean(mean(full(S(S < 0 & non_diag_mask)))));
            fprintf("Positive S min-max: [%.5f, %.5f]\n", ...
                    min(min(full(S(S > 0 & non_diag_mask)))), max(max(full(S(S > 0 & non_diag_mask)))));
        end

        if threshold == 0
            % use the mean of non-diagonal (hence positive) values
            threshold = mean(S(S>0&non_diag_mask));
        end
        
        % S has negative values only on diagonal, hence the
        %   threshold is just for positive values
        S(non_diag_mask & S < threshold) = 0;
        
        if verbose
            fprintf("post mask density(S): %2.5f\n", (nnz(S)*100)/numel(S));
        end
    end
    
    % create D_s and C
    D_s = sqrt(D);
    C = chol(-S);
    
    if verbose
        fprintf("density(C): %2.5f\n", (nnz(C)*100)/numel(C));
    end
    
    % build the full preconditioner matrix if needed
    if build_full_precon
        Z = sparse(size(C, 1), n_edges);
        P = [D Z'; Z -S];
    end
end
