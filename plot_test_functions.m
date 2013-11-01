function [ ] = plot_test_functions(P, xn, qd, np)
  % function [] = plot_test_functions(Vopt, qd, np)
  %
  % PURPOSE: Plots each component of the optimal test functions on separate
  % plots. Vopt contains a different test function in each of its pu+pq+2
  % columns. The upper pw+1 rows of each column contain the first state 
  % component of the test function, while the lower pv+1 rows containt the 
  % second component.
  %
  % INPUTS:
  %   P : projection to optimal lower-order test functions
  %   xn : basis nodes
  %   qd : quadrature data
  %   np : number of plotting points
  %
  % OUTPUTS:
  %

  % test function basis nodes
  nn = size(P, 2);
  p = nn-1;

  % reference space for plotting
  xi = linspace(0, 1, np)';

  % interpolate test functions to plotting points
  Phi = basis(xn, xi);
  v = Phi*P;

  % plot each component of all test functions
  plot(xi, v, 'LineWidth',2); xlabel('\xi'); ylabel('Test functions');

end

