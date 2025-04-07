function h=pltellips(mu,Q,varargin)
% pltellips   plot an error ellipse
%    PLTELLIPS(MU,Q) plots the error ellips associated to the 2x2 covariance
%    matrix Q at position MU.
%
%    PLTELLIPS(MU,Q,...) other arguments will be passed to the plot function.
%
%    H=PLTELLIPS(...) returns the plot handle H.
%
%    See also plot. 

if nargin < 2 
   error('not the correct number of input arguments')
end

% compute eigenvalues and eigen vectors

[v,d] = eig(Q); % Compute eigen-stuff

% Check that the matrix is positive definite
if any(diag(d) <=0)
  error('The covariance matrix must be positive definite')
end

% generate an eliips with n points

n=100;
p=0:pi/n:2*pi;
ell = [cos(p'),sin(p')] * sqrt(d) * v'; 

hold_state = get(gca,'nextplot');
h=plot(mu(1)+ell(:,1),mu(2)+ell(:,2),varargin{:});
set(gca,'nextplot',hold_state);

end