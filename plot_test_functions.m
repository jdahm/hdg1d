function [ ] = plot_test_functions(Vopt, qd)

% function [] = plot_test_functions(Vopt, qd)
%
% PURPOSE: Plots each component of the optimal test functions on separate
% plots. Vopt contains a different test function in each of its pu+pq+2
% columns. The upper pw+1 rows of each column contain the first state 
% component of the test function, while the lower pv+1 rows containt the 
% second component.
%
% INPUTS:
%   Vopt: optimal test function coefficients ([pv+pw+2]x[pu+pq+2])
%     qd: quadrature data
%
% OUTPUTS: None

% test function basis nodes
xnv = qd.xnv;
xnw = qd.xnw;
pw = qd.pw;

% reference space for plotting
Xiplot = [0:.005:1]'; 

% interpolate test functions to plotting points
vPhiPlot = basis(xnv, Xiplot);
wPhiPlot = basis(xnw, Xiplot);
wPlot = wPhiPlot*Vopt(1:pw+1,:);
vPlot = vPhiPlot*Vopt(pw+2:end,:);

% plot each component of all test functions
figure;
plot(Xiplot, wPlot,'LineWidth',2); xlabel('Xi'); ylabel('Test functions'); title('w component');
hold off; figure;
plot(Xiplot, vPlot,'LineWidth',2); xlabel('Xi'); ylabel('Test functions'); title('v component');

end

