problem_size = 500;
A = generate_graph_matrix(problem_size, 0.5);
b = rand(problem_size, 1);
t=1:1:problem_size;
num_runs = 10;
qr_tot_new_times = zeros(1, problem_size);
qr_tot_old_times = zeros(1, problem_size);
mat_mul2_tot_new_times = zeros(problem_size, 1);
mat_mul2_tot_old_times = zeros(problem_size, 1);

for run_index = 1:num_runs
    fprintf("Run number %g\n", run_index);
    [x_new, res_new, setup_time, lanczos_times, qr_times_new, inv_times_opt, final_times, mat_mul1_times, mat_mul2_times_new] = minres_qr(A, b, true);
    [x_old, res_old, setup_time, lanczos_times, qr_times_old, inv_times_old, final_times, mat_mul1_times, mat_mul2_times_old] = minres_qr(A, b, false);
    qr_tot_new_times = qr_tot_new_times + qr_times_new;
    qr_tot_old_times = qr_tot_old_times + qr_times_old;
    mat_mul2_tot_new_times = mat_mul2_tot_new_times + mat_mul2_times_new;
    mat_mul2_tot_old_times = mat_mul2_tot_old_times + mat_mul2_times_old;
end


qr_tot_new_times = qr_tot_new_times / num_runs;
qr_tot_old_times = qr_tot_old_times / num_runs;
mat_mul2_tot_new_times = mat_mul2_tot_new_times / num_runs;
mat_mul2_tot_old_times = mat_mul2_tot_old_times / num_runs;
size(qr_tot_new_times)
size(mat_mul2_tot_new_times)

x_diff = norm(x_new - x_old);
fprintf("The difference between the solutions is %d\n", x_diff);

subplot(3, 2, 1);
plot(t, mat_mul2_tot_new_times, t, qr_tot_new_times);
legend("mul2", "total");
title("new");

subplot(3, 2, 2);
plot(t, mat_mul2_tot_old_times, t, qr_tot_old_times);
%plot(t, log(mat_mul2_times_old), t, log(qr_times_old));
legend("mul2", "total");
title("old");

subplot(3, 2, 3);
plot(t, qr_tot_old_times, t, qr_tot_new_times);
%plot(t, log(mat_mul2_times_old), t, log(qr_times_old));
legend("old", "new");
title("qr comparison");

subplot(3, 2, 4);
plot(t, mat_mul2_tot_old_times, t, mat_mul2_tot_new_times);
%plot(t, log(mat_mul2_times_old), t, log(qr_times_old));
legend("old", "new");
title("mul comparison");

subplot(3, 2, 5);
mat_mul_diff = mat_mul2_tot_new_times - mat_mul2_tot_old_times;
plot(t, mat_mul_diff);
%plot(t, log(mat_mul2_times_old), t, log(qr_times_old));
%legend("mul2");
title("diff");

subplot(3, 2, 6);
plot(t, inv_times_old, t, inv_times_opt);
%plot(t, log(mat_mul2_times_old), t, log(qr_times_old));
legend("old", "new");
title("inv comparison");

old_qr_mean = mean(qr_tot_old_times);
new_qr_mean = mean(qr_tot_new_times);
fprintf("Old QR mean: %d\n", old_qr_mean);
fprintf("New QR mean: %d\n", new_qr_mean);

old_mul_mean = mean(mat_mul2_tot_new_times(2:end));
new_mul_mean = mean(mat_mul2_tot_old_times(2:end));
fprintf("Old mul mean: %d\n", old_mul_mean);
fprintf("New mul mean: %d\n", new_mul_mean);
