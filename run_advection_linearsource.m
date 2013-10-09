% Advection-only case with linear source term

fd.a = 1e-1;
fd.b = 0;
fd.c = 1.0;
fd.q_present = false;
fd.stab_type = 'upwind';
fd.vl = 1;

sd.present = true;
sd.type = 'linear';
sd.c = 1.0;

%dt = 1.0;
%td.c = 1.0/dt;
td.c = 0.0;
td.u_t = 0.0;

pq = 4;
pv = 4;
pu = 4;
pw = 4;

md = mesh(0., 2., 3);

xnq = create_nodes(pq, 'SegLagrangeGauss');
xnv = create_nodes(pv, 'SegLagrangeGauss');
xnu = create_nodes(pu, 'SegLagrangeGauss');
xnw = create_nodes(pw, 'SegLagrangeGauss');

lbd.type = 'd';
lbd.data = 1;
rbd.type = 'n';
rbd.data = 0;

[Q, U, L] = initialize(pq, pu, md.ne);

qd = quadrature_data(xnq, xnu, xnv, xnw);

[dQ, dU, dL] = hdg_solve(Q, U, L, lbd, rbd, md, td, fd, sd, qd);

Q = Q + dQ;
U = U + dU;
L = L + dL;

figure(1); clf;
h1 = plot_elems(md.xs, md.xe, xnu, U, 100);

% plot exact solution
x = linspace(0,1,100);
y = exp(-sd.c/fd.a*x);
hold on;
plot(x,y,'--');
hold off;
