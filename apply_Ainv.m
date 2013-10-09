function [Q, U] = apply_Ainv(AQQ, AQU, AUQ, AUU, Q, U, q_present)
  % function [Q, U] = apply_Ainv(AQQ, AQU, AUQ, AUU, Q, U, q_present)
  %
  % PURPOSE: Takes a pre-processed block matrix A and performs
  %   A^{-1}*[Q;U]
  % where Q and U can either form a matrix or vector.
  %
  % INPUTS:
  %   {AQQ, AQU, AUQ, AUU} : matrix blocks after pre-processing
  %   Q : block [nq, c]
  %   U : block [nu, c]
  %   q_present : indicator of number of blocks
  %
  % OUTPUTS:
  %   Q : first nq rows of A^{-1}*[Q;U]
  %   U : last nu rows of A^{-1}*[Q;U]
  %

  if q_present
    % needed since AQQ should not be invertible unless q is present
    % Q = AQQ^{-1}*Q
    Q = AQQ\Q;
  end
  % W = -AUU^{-1}*U
  W = -AUU\U;
  if q_present
    % W += AUQ*Q
    W = W+AUQ*Q;
    % Q += AQU*W
    Q = Q+AQU*W;
  else
    Q = zeros(size(Q));
  end
  % U = -W
  U = -W;
