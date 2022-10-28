function [b] = multiply_preconditioner(a, D_s, C, transpose_C)
% Indirectly multiplies the Schur Complement preconditioner P.
%   What is being computed is the following: given the factor B such that
%   B'B = chol(P) and B = [D_s 0; 0 C], b = B^{-1} a. 
%   This is achieved computing the above multiplication via individual blocks.
%   Considering a = [a_1; a_2] and b = [b_1; b_2], 
%   b_1 = D_s^{-1} * a_1 and b_2 = C^{-1} * a_2.
%   
%
% Inputs:
% 
% - a ([nodes+edges x 1[ real vector): the vector that has to be multiplied
%       by the preconditioner.
% 
% - D_s ([edges x edges] real matrix): the square root of the D.
%
% - C ([nodes x nodes] real matrix): the Cholesky factor of -S.
%
% Optional inputs:
%
% - transpose_C (bool, defaults to false): whether to transpose C or not, 
%       thus transposing the whole preconditioner, as D_s is diagonal and 
%       the other two blocks are 0.
%       This is used just for clarity, the same effect can obviously be
%       achieved by passing C'.
%
% Outputs:
% 
% - b ([nodes+edges x 1[ real vector): the preconditioned vector b = B^{-1} a.
%
    % handle optional parameter
    if exist('transpose_C', 'var') == 0
       transpose_C = false;
    end
    
    % compute b_1 by solving the diag system D_s*{-1} * a_1
    size_D_s = size(D_s, 1);
    b_1 = diag_inverter(D_s, a(1:size_D_s, :)); % D_s
    
    % compute b_2 by solving C^{-1} * a_2
    if transpose_C == true
       C = C';
    end
    b_2 = C \ a(size_D_s+1:end, :);
    
    b = [b_1; b_2];
end

function [x] = diag_inverter(D, b)
% Directly solves the diagonal system D_s x = b by computing the inverse of
%   each diagonal element and multiplying it element-wise to b. This is
%   faster than Matlab mldivide (\) operator, as it does not have a case
%   for diagonal matrices.
    if issparse(D)
        x = D \ b;
    else
        x = 1./diag(D) .* b;
    end
end
