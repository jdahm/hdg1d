% Quadratic manufactured solution

scheme = 'hdg';
%scheme = 'dpg';

fd.a = 1;
fd.b = 1e-9;
fd.q_present = true;
fd.stab_type = 'upwind';
fd.vl = 1.0;
fd.c = 1.0;

sd.present = true;
sd.type = 'geometric';
sd.name = 'ms_sine';

td.c = 0.0;
td.u_t = 0.0;

pq = 3;
pu = 3;

pv = 3;
pw = 3;

% optimal test function weights (higher weight means more emphasis on
% boundary portion of output)
testd.wleft = 1e6;
testd.wright = 1e6;

md = mesh(0., 1., 32);

qbasis = 'SegLagrange';
ubasis = 'SegLagrange';

xnq = create_nodes(pq, qbasis);
xnv = create_nodes(pv, qbasis);
xnu = create_nodes(pu, ubasis);
xnw = create_nodes(pw, ubasis);

lbd.type = 'd';
lbd.data = 0.;
rbd.type = 'd';
rbd.data = 0.;

[Q, U, L] = initialize(pq, pu, md.ne);

qd = quad_data(xnq, xnu, xnv, xnw);

if strcmp(scheme, 'dpg')

  % compute optimal test functions
  Vopt = compute_optimal_test_functions(fd, qd, md, sd, lbd, rbd, testd);

  % plot test functions
  figure(3); clf;
  plot_test_functions(Vopt.w, xnw, qd, 1000);
  figure(4); clf;
  plot_test_functions(Vopt.v, xnv, qd, 1000);

  % fill DPG quadrature data
  qd = dpg_quad_data(qd, Vopt.v, Vopt.w);
end

[dQ, dU, dL] = solve_linear_system(scheme, Q, U, L, lbd, rbd, ...
				   md, td, fd, sd, qd);

Q = Q + dQ;
U = U + dU;
L = L + dL;

figure(1); clf;
h1 = plot_elems(md.xs, md.xe, xnu, U, 1000);
if md.ne > 1
  h1 = plot_traces(md.xs, md.xe, L);
end

figure(2); clf;
h1 = plot_elems(md.xs, md.xe, xnq, Q, 1000);

% exact solution on top
figure(1);
x = linspace(md.xs, md.xe, 1000);
y = sin(2*pi*x);
hold on;
plot(x,y,'--');
hold off;

figure(2);
y_x = 2*pi*cos(2*pi*x);
hold on;
plot(x,y_x,'--');

%errl = error_norm(Q, U, L, xnq, xnu, pu+1, md, fd, ...
%		  'Trace', @(fd, x) sin(2*pi*x))
%erru = error_norm(Q, U, L, xnq, xnu, pu+1, md, fd, ...
%		  'L2', @(fd, x) sin(2*pi*x))
%errq = error_norm(Q, U, L, xnq, xnu, pu+1, md, fd, ...
%		  'H1', @(fd, x) sin(2*pi*x), @(fd, x) 2*pi*cos(2*pi*x))
