% This function returns the steady-state translation rate given \lambda_0,...\lambda_n.
% It constructs the matrix A defined in the RFM concavity paper, computes its eigenvalues and 
% returns R = 1/(max(eigenvalues))^2.
%
%   Usage: [ R, e, A ] = RFM_n_R_eval ( l )
%  
%   Where:  l - an n+1 size vector of rates [\lambda_0,...,\lambda_n]
%
%           e - a vector of A's eigenvalues
%           R - steady-state translation rate
%           A - the matrix defined in the RFM concavity paper
%
%
% Yoram Zarai, 5/29/14

% ======================================================================================================

function [ R, e, A ] = RFM_n_R_eval ( l )

  % if one (or more) l's is zero, then R=0
  if( sum( l == 0 ) )
    e = [];
    A = [];
    R = 0;
    return;
  end;
  
q = l .^ (-0.5); % \lambda_i^(-0.5)
A = diag( q, -1 ) + diag( q, 1 );
e = eig( A );
R = max( e ) ^(-2); % A is symmetric, thus eigenvalues are real

% example for n=4
%
% A = [ 0      q(1)      0      0      0      0 ;...
%       q(1)   0         q(2)   0      0      0 ;...
%       0      q(2)      0      q(3)   0      0 ;...
%       0      0         q(3)   0      q(4)   0 ;...
%       0      0         0      q(4)   0      q(5) ;...
%       0      0         0      0      q(5)   0 ];

