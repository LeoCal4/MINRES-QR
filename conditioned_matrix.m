function A=conditioned_matrix(m,epsilon)

% Generates a nonsingular matrix while bounding the condition number
%
% Call sequence: A = TestMatrix(m,epsilon)
%
% INPUT:
%    m        the dimension of the matrix
%    epsilon  a real number in (0,1)
%
% OUTPUT:
%    A        an m by m nonsingular matrix
%
% The infinity norm condition number of A is bounded by 
% 
%                (1+epsilon)/(1-epsilon)
% 

% Programming by Carl Christian Kjelgaard Mikkelsen
%   2021-10-22 Initial programming and testing

% Generate a random matrix with entries uniformly distributed in (0,1)
E=rand(m,m);
% Compute the sum of the absolute values of E along each row
w=abs(E)*ones(m,1);
% Scale E to ensure that the infinity norm is epsilon
E=epsilon*(diag(w)\E);
% Construct the final matrix
A=eye(m,m)-E;