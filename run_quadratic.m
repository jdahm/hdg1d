% Quadratic manufactured solution

fd.a = 1;
fd.b = 1;
fd.q_present = true;
fd.stab_type = 'upwind';
fd.vl = 1.0;

sd.present = true;
sd.type = 'ms_quadratic';

%dt = 1.0;
%td.c = 1.0/dt;
td.c = 0.0;
td.u_t = 0.0;

pq = 2;
pv = 2;
pu = 2;
pw = 2;

md = mesh(0., 2., 3);

xnq = create_nodes(pq, 'SegLagrange');
xnv = create_nodes(pv, 'SegLagrange');
xnu = create_nodes(pu, 'SegLagrange');
xnw = create_nodes(pw, 'SegLagrange');

lbd.type = 'd';
lbd.data = 0.;
rbd.type = 'd';
rbd.data = 0.;

[Q, U, L] = initialize(pq, pu, md.ne);

qd = quadrature_data(xnq, xnu, xnv, xnw);

[dQ, dU, dL] = hdg_solve(Q, U, L, lbd, rbd, md, td, fd, sd, qd);

Q = Q + dQ;
U = U + dU;
L = L + dL;

figure(1);
h1 = plot_elems(md.xs, md.xe, xnu, U, 1000);
h1 = plot_traces(md.xs, md.xe, L);

% exact solution
x = linspace(0, 2, 1000);
y = -x.*(x-2);
hold on;
plot(x,y,'--');
hold off;
