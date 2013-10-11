function xn = create_nodes(p, basis_type)

  switch (basis_type)
    case 'SegLagrange'
      xn = linspace(0, 1, p+1)';
    case 'SegLagrangeGauss'
      [xn,~] = lgwt(p+1,0,1);
    otherwise
      error('unknown basis type');
  end
