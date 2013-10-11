% Boundary layer solution

fd.a = 1;
fd.b = 1e-2;
fd.q_present = true;
fd.stab_type = 'upwind';
fd.vl = fd.a/fd.b;
%fd.vl = 1.;

sd.present = false;

%dt = 1.0;
%td.c = 1.0/dt;
td.c = 0.0;
td.u_t = 0.0;

pq = 5;
pv = 5;
pu = 5;
pw = 5;

md = mesh(0., 1., 2);

xnq = create_nodes(pq, 'SegLagrange');
xnv = create_nodes(pv, 'SegLagrange');
xnu = create_nodes(pu, 'SegLagrange');
xnw = create_nodes(pw, 'SegLagrange');

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

figure(2); clf;
h1 = plot_elems(md.xs, md.xe, xnu, U, 100);
h1 = plot_traces(md.xs, md.xe, L);

% plot exact solution on top
x = linspace(0,1,1000);
c1 = fd.a/fd.b/(exp(fd.a/fd.b)-1);
c2 = -1/(exp(fd.a/fd.b)-1);
y = fd.b/fd.a*c1*exp(fd.a/fd.b*x)+c2;

hold on;
plot(x, y, '--');
