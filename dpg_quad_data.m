function dqd = dpg_quad_data(qd, Mv, Mw)
  % function dqd = dpg_quad_data(qd, Mv, Mw)
  %
  % PURPOSE: computes a new quad data struct with "optimal" test
  % functions.
  %
  % INPUTS:
  %   qd : quad data with regular polynomial test functions [struct]
  %   Mv : matrix of optimal test function coefficients for v
  %   Mw : matrix of optimal test function coefficients for w
  %
  % OUTPUTS:
  %   dqd : quad data struct with optimal test functions [struct]
  %

  % copy quad data struct
  dqd = qd;

  % test bases
  dqd.vPhi = qd.vPhi * Mv;
  dqd.vPhi0 = qd.vPhi0 * Mv;
  dqd.vPhi1 = qd.vPhi1 * Mv;

  dqd.wPhi = qd.wPhi * Mw;
  dqd.wPhi0 = qd.wPhi0 * Mw;
  dqd.wPhi1 = qd.wPhi1 * Mw;

  % reference-space test function gradients
  dqd.vGPhi = qd.vGPhi * Mv;
  dqd.wGPhi = qd.wGPhi * Mw;
