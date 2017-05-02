% This function implements the Ordinary Differential Equations of RFM
%
% The following global variables must be defined:
% 1. RFM_n - the value of n (dimension of RFM). Must be >= 1.
% 2. RFM_lm - an n+1 vector containing the \lambdas:  
%             RFM_lm = [ \lambda_0, \lambda_1,...,\lambda_n ];
%
%
% Yoram, Zarai, 12/8/13

% ================================================================================

function ret = rfm_ode( t, x_vec )

  global RFM_n RFM_lm
  
  % using x, all ODE equations have the same structure
  x = [ 1 ; x_vec(:) ; 0 ]; % so x(0)=1, x(n+1)=0, x(2:end-1) = x_vec
  ret = zeros( RFM_n, 1 );
  for i = 1 : RFM_n
    ret( i ) = RFM_lm(i) * x(i) * ( 1 - x(i+1) ) - RFM_lm(i+1) * x(i+1) * ( 1 - x(i+2) );
  end;
  
  