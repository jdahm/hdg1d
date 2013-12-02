function err = error_norm(Q, U, L, xnq, xnu, qorder, md, fd, sd, norm_type, varargin)

  % perform a sanity check and pull off the anonymous functions in varargin
  switch norm_type
    case {'L1', 'L2', 'Trace', 'CellAvgU'}
      if length(varargin) < 1
	error('need u(fd,sd,x)');
      end
      u = varargin{1};
    case 'CellAvgQ'
      if length(varargin) < 1
	error('need u_x(fd,sd,x)');
      end
      u_x = varargin{1};
    case 'H1'
      if length(varargin) < 2
	error('need u(fd,sd,x), u_x(fd,sd,x)');
      end
      u = varargin{1};
      u_x = varargin{2};
    otherwise
      error('unknown norm type');
  end

  dx = (md.xe - md.xs) / md.ne;

  [x,w] = lgwt(qorder, 0, 1);
  qPhi = basis(xnq, x);
  uPhi = basis(xnu, x);

  nnq = length(xnq);
  nnu = length(xnu);

  % loop over elements
  err = 0.;
  erru = 0.;
  errq = 0.;
  for elem=1:md.ne

    % global quad point locations
    xg = md.xs+(elem-1)*dx+dx*x;

    % get the approximate solution at the quad points
    if strcmp(norm_type, 'Trace')
       uhle = u(fd, (elem-1)*md.dx);
       uhre = u(fd, elem*md.dx);
       if elem ~= 1
	 uhla = L(elem-1);
       else
	 uhla = uhle;
       end
       if elem ~= md.ne
	 uhra = L(elem);
       else
	 uhra = uhre;
       end
    end
    if strcmp(norm_type, 'L1') || strcmp(norm_type, 'L2') || strcmp(norm_type, 'H1') || strcmp(norm_type, 'CellAvgU')
      UE = U((elem-1)*nnu+1:elem*nnu);
      ua = uPhi*UE;
      ue = u(fd, sd, xg);
    end
    if strcmp(norm_type, 'H1') || strcmp(norm_type, 'CellAvgQ')
       QE = Q((elem-1)*nnq+1:elem*nnq);
       u_xa = qPhi*QE;
       u_xe = u_x(fd, sd, xg);
    end

    % add to the norms
    switch norm_type
      case 'Trace'
	err = err + md.dx*((uhla - uhle)^2 + (uhra - uhre)^2);
      case 'L1'
	err = err + dx*sum(w.*abs(ua-ue));
      case 'L2'
	err = err + dx*sum(w.*(ua-ue).^2);
      case 'H1'
	erru = erru + dx*sum(w.*(ua-ue).^2);
	errq = errq + dx*sum(w.*(u_xa-u_xe).^2);
      case 'CellAvgU'
	err = err + dx*abs(sum(w.*ua) - sum(w.*ue));
      case 'CellAvgQ'
	err = err + dx*abs(sum(w.*u_xa) - sum(w.*u_xe));
      otherwise
	error('unknown norm type');
    end

  end

  % post-summation operation
  if strcmp(norm_type, 'L2') || strcmp(norm_type, 'Trace')
    err = sqrt(err);
  elseif strcmp(norm_type, 'H1')
    err = sqrt(erru + errq);
  end
