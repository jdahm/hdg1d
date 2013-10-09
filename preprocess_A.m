function [AQQ, AQU, AUQ, AUU] = preprocess_A(AQQ, AQU, AUQ, AUU, q_present)
  % function [AQQ, AQU, AUQ, AUU] = preprocess_A(AQQ, AQU, AUQ, AUU, q_present)
  %
  % PURPOSE: Takes a block matrix A
  %   A = [AQQ, AQU; AUQ, AUU]  if q_present is true
  %   A = [AUU]                 otherwise
  % and prepares it to be quickly inverted.
  %
  % INPUTS:
  %   {AQQ, AQU, AUQ, AUU} : matrix blocks
  %   q_present : indicator of number of blocks
  %
  % OUTPUTS:
  %   AQQ : stays same (AQQ)
  %   AQU : becomes AQQ^{-1}*AUQ
  %   AUQ : K^{-1}*AUQ where K=(AUU-AUQ*AQQ^{-1}*AUQ)
  %   AUU : K^{-1}
  %

  if q_present % needed because AQQ should not be invertible unless q is present
    % AQU = AQQ^{-1}*AQU
    AQU = AQQ\AQU;
    % AUU -= AUQ*AQU
    AUU = AUU-AUQ*AQU;
    % AUQ = AUU^{-1}*AUQ
    AUQ = AUU\AUQ;
  end
