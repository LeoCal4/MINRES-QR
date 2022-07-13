A = rand(2, 1);
prev_Q = eye(1);
for iter = 1:500
    [prev_Q, R_k] = iterative_QR(A, prev_Q, iter);
    [Q, R] = qr(A);
    if norm(Q-prev_Q) > 1e-13 || norm(R-R_k) > 1e-13
       disp("cazzo")
       norm(Q-prev_Q)
       norm(R-R_k)
       iter
       return
    end
    A = expand_matrix_as_lanczos(A);
end
disp("ok :)")
norm(Q-prev_Q)
norm(R-R_k)