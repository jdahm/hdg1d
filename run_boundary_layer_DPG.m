% Boundary layer solution
close all; clc; clear all;

% flux data 
fd.a = 1.2;
fd.b = .1;
fd.q_present = true;
fd.stab_type = 'upwind';
fd.vl = 1e5;

% source data
sd.present = true;
sd.type = 'linear';
sd.a = 0;
sd.b = 0;

%dt = 1.0;
%td.c = 1.0/dt;
td.c = 0.0; % time data (not used for now)
td.u_t = 0.0;

% orders of trial and test spaces (both u and q)
pu = 2;
pq = 2;
pw = 7;
pv = 7;

% optimal test function weights (higher weight means more emphasis on
% boundary portion of output)
testd.wleft = 1000;
testd.wright = 1000;

% mesh extent and number of elements
md = mesh(0., 1., 1);

% basis node locations corresponding to specified orders
xnq = create_nodes(pq, 'SegLagrange');
xnv = create_nodes(pv, 'SegLagrange');
xnu = create_nodes(pu, 'SegLagrange');
xnw = create_nodes(pw, 'SegLagrange');

% boundary conditions
lbd.type = 'd'; % left boundary
lbd.data = 0.;
rbd.type = 'd'; % right boundary
rbd.data = 1;

% initialize solution to zeros
[Q, U, L] = initialize(pq, pu, md.ne);

% fill quadrature data
qd = quad_data(xnq, xnu, xnv, xnw);

% compute optimal test functions
Vopt = compute_optimal_test_functions(fd, qd, md, sd, testd);

% plot test functions
plot_test_functions(Vopt.wv, qd);

% fill DPG quadrature data
qd = dpg_quad_data(qd, Vopt.v, Vopt.w);

% obtain updates to q, u, and trace variables
%[dQ, dU, dL] = hdg_solve(Q, U, L, lbd, rbd, md, td, fd, sd, qd);

% do HDPG solve
[dQ, dU, dL] = hdgdpg_solve(Q, U, L, lbd, rbd, md, td, fd, sd, qd);

% update solutions
Q = Q + dQ;
U = U + dU;
L = L + dL;

figure;
h1 = plot_elems(md.xs, md.xe, xnu, U, 100);
%h1 = plot_traces(md.xs, md.xe, L);

% exact solution
if fd.a ~=0
    y = @(fd, x) (exp(fd.a/fd.b*x) - 1.0)./(exp(fd.a/fd.b)-1);
    y_x = @(fd, x) fd.a/fd.b*exp(fd.a/fd.b*x)./(exp(fd.a/fd.b)-1);
else
    y = @(fd, x) lbd.data*(1 - x) + x*rbd.data;
    y_x = @(fd, x) rbd.data - lbd.data;
end

% plot exact solution on top
x = linspace(0,1,1000);
hold on;
plot(x, y(fd, x), '--');
title('solution, u')
hold off;

% plot numerical gradient, q
figure;
h2 = plot_elems(md.xs, md.xe, xnq, Q, 100);

% plot exact gradient
hold on;
plot(x, y_x(fd, x), 'r--');
title('gradient, q');
hold off;

% get error norm (L1, L2, or H1)
err = error_norm(xnq, xnu, Q, U, L, 10, md, fd, 'L2', y)
