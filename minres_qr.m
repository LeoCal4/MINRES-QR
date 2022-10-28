function [x, residuals, iter] = minres_qr(A_t, b_t, precon, n_edges, precon_threshold, reorthogonalize, verbose)
% Computes the full MINRES-QR algorithm, starting from sparse matrix A_t and
%   vector b_t.
% At each iteration, the solution to the system is searched in an
%   increasingly growing Lanczos subspace, whose bases are iteratively
%   obtained by the Lanczos algorithm. The tridiagonal matrix T obtained
%   at iteration k is then decomposed using an iterative version of QR
%   factorization, which updates the matrices Q and R at each step,
%   instead of recomputing them from scratch. The final solution for each 
%   iteration is finally obtained by solving the system R_k x_k = c, where
%   c is a vector obtained from of \beta_k Q_k.
% The algorithm stops once the relative residual reaches the desired
%   precision (1e-09).
% Additional parameters are available to precondition and
%   reorthogonalize the problem at hand.
% In case of (Schur Complement) preconditioning, the preconditioner factors D_s and C are
%   obtained from A_t; they are then applied to b_t, in order to obtain its
%   preconditioned version. At each iteration, they are passed to the
%   Lanczos function, which accordingly preconditions the matrix A_t when it
%   is applied to v_k.
%
% Inputs:
%
% - A_t ([(nodes-1)+edges x (nodes-1)+edges] real sparse matrix): the truncated
%       version of the problem matrix. It is composed by [D F'; F 0], where D is a diagonal
%       positive definite matrix, while F is the incidence matrix of the graph 
%       with one row removed.
%
% - b_t ([(nodes-1)+edges x 1] real column vector): the truncated version of
%       the original b vector
%
% Optional inputs:
%
% - precon (bool): if true, the system is preconditioned using Schur
%       Complement preconditioner.
%
% - n_edges (int, scalar): the number of edges, needed to extract the D
%       matrix from A_t during the creation of the preconditioner.
%
% - precon_threshold (real, scalar, defaults to -1): specifies the desired 
%       non-negative threshold for non-diagonal elements in the S block.
%       It can be defined in the two following ways:
%       - threshold > 0: the desired threshold value;
%       - threshold = 0: the threshold value is set to be the mean of the
%           non-diagonal values of S.
%
% - reorthogonalize (bool, defaults to false): whether to reorthogonalize
%       the system or not during Lanczos step. This noticebly increases its cost,
%       while making the algorithm backward stable.
%
% - verbose (bool, defaults to false): whether to print information about
%       the algorithm or not.
%
% Outputs:
%
% - x ([(nodes-1)+edges x 1] real vector): the computed solution of the system A_t x_t = b_t.
%       In order to obtain the solution to the real system Ax = b, simply
%       append a 0 at the end of the vector.
%
% - residuals ([iter x 1] real vector): the relative residuals of all the
%       iterations.
%
% - iter (int, scalar): the number of iterations needed for the algorithm
%       to converge.
%

    % handle optional variables
    if exist('precon', 'var') == 0
       precon = false;
    end
    if precon && exist('n_edges', 'var') == 0
        error("ERROR: must specify the number of edges in the graph/the dimension of the D matrix in order to apply the Schur Preconditioner");
    end
    if exist('precon_threshold', 'var') == 0
        precon_threshold = -1;
    end
    if exist('reorthogonalize', 'var') == 0
       reorthogonalize = false;
    end
    if exist('verbose', 'var') == 0
       verbose = false;
    end
    
    % saving the original b_t, before eventual preconditioning
    b_start = b_t;
    
    % initializing preconditioning matrices to default values
    D_s = 0;
    C = 0;
    if precon == true
        if verbose; fprintf("\tInitializing preconditioner\n"); end
        [D_s, C] = create_preconditioner(A_t, n_edges, precon_threshold, false, verbose);
        b_t = multiply_preconditioner(b_t, D_s, C, true);
    end
    
    if reorthogonalize && verbose
        fprintf("\tApplying reorthogonalization during Lanczos\n");
    end
    
    % initialized general vars
    size_A = size(A_t, 1);
    threshold = 1e-9;
    residuals = nan(1, size_A);
    
    % initialize Lanczos' matrices
    %   preinitializing them with the full dimensions makes them too big
    %   for problems with notable sizes (changes in times are minimal)
    V = zeros(size_A, 2);
    T = zeros(2, 1);
    prev_w = b_t;
    beta_1 = norm(b_t);
    
    % initialize QR and final solution matrices
    Q = eye(1);
    R = eye(1);
    c = zeros(1, 1);
    
    if verbose; fprintf("Starting iterations\n"); end
    
    % iterate up to size of A_t at max
    for k = 1:size_A
        if verbose; fprintf("Iteration %i\n", k); end
        
        % LANCZOS
        [V, T, prev_w] = iterative_lanczos(A_t, V, T, prev_w, k, reorthogonalize, precon, D_s, C);
        
        % QR
        [Q, R] = iterative_QR(T(1:k+1, 1:k), Q, R, k);
        
        % BACKSUBSTITUTION
        %  use the reduced version of the QR, using just Q_1 and R_1
        %      from Q = [Q_1 Q_2] and R = [R_1; 0]
        c(k, 1) = beta_1 * Q(1, k);
        y_k = R(1:end-1, :) \ c;
        
        % SOLUTION
        x_k = V(:, 1:k) * y_k;
        if precon == true
            % solution of the original truncated system
            x_k = multiply_preconditioner(x_k, D_s, C); % x = B^-1 * \hat(x)
        end
        
        % save the current relative residual
        residual = norm(b_start - A_t * x_k) / norm(b_start);
        residuals(k) = residual;
        
        % check whether the residual is small enough
        if residual <= threshold
            break
        end
    end
    
    x = x_k;
    residuals = residuals(1:k);
    if verbose; fprintf("MINRES-QR converged in %i iterations\n", k); end
    iter = k;
end
