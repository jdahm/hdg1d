% Boundary layer solution

fd.a = 1;
fd.b = 1e-1;
fd.q_present = true;
fd.stab_type = 'upwind';
fd.vl = fd.a/fd.b;
%fd.vl = 1.;

sd.present = false;

%dt = 1.0;
%td.c = 1.0/dt;
td.c = 0.0;
td.u_t = 0.0;

pq = 2;
pv = 2;
pu = 2;
pw = 2;

md = mesh(0., 1., 4);

xnq = create_nodes(pq, 'SegLagrangeGauss');
xnv = create_nodes(pv, 'SegLagrangeGauss');
xnu = create_nodes(pu, 'SegLagrangeGauss');
xnw = create_nodes(pw, 'SegLagrangeGauss');

lbd.type = 'd';
lbd.data = 0.;
rbd.type = 'd';
rbd.data = 1.;
 
[Q, U, L] = initialize(pq, pu, md.ne);

qd = quadrature_data(xnq, xnu, xnv, xnw);

[dQ, dU, dL] = hdg_solve(Q, U, L, lbd, rbd, md, td, fd, sd, qd);

Q = Q + dQ;
U = U + dU;
L = L + dL;

figure(1); clf;
h1 = plot_elems(md.xs, md.xe, xnu, U, 100);
%h1 = plot_traces(md.xs, md.xe, L);

% exact solution
y = @(fd, x) (exp(fd.a/fd.b*x) - 1.0)./(exp(fd.a/fd.b)-1);
y_x = @(fd, x) fd.a/fd.b*exp(fd.a/fd.b*x)./(exp(fd.a/fd.b)-1);

% plot exact solution on top
x = linspace(0,1,1000);
hold on;
plot(x, y(fd, x), '--');
hold off;

% get error norm
err = error_norm(xnq, xnu, Q, U, L, 10, md, fd, 'H1', y, y_x)
