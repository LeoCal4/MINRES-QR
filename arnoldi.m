function [Q, H]  = arnoldi(A, b, n)
% returns a matrix Q of basis of K_{n+1}(A,b)

m = length(b);
Q = zeros(m, n+1);
Q(:, 1) = b / norm(b);
H = zeros(n+1, n);
for j = 1:n
   w = A*Q(:,j);
   for i = 1:j
       H(i, j) = Q(:,i)'*w;
      w = w - Q(:,i)*H(i, j);
   end
   H(j+1, j) = norm(w);
   Q(:,j+1) = w / H(j+1, j);
end