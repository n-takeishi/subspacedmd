function [lam, w, z] = orddmd(Y0, Y1, r)
% ORDDMD Compute (exact) dynamic mode decomposition based on SVD.
%
% <Input>
%   Y0, Y1: N x M snapshot matrices, where N is the dimensionality
%     of snapshots and M is the number of snapshots.
%   r: Number of retained dynamic modes. If omitted, determined automatically.
%
% <Output>
%   lam: Eigenvalues of matrix A, where A is an approximation of Koopman
%     operator.
%   w  : Right eigenvectors of matrix A.
%   z  : Left eigenvectors of matrix A. The values of eigenfunctions of
%     Koopman operator is computed by z'*Y.
%
% For details, see the following paper:
%   Schmid, Dynamic mode decomposition of numerical and experimental data,
%     Journal of Fluid Mechanics, 656:5-28, 2010.
%   Tu et al., On dynamic mode decomposition: theory and applications,
%     Journal of Computational Dynamics, 1(2):391-421, 2014.

if nargin<3, r=rank(Y0); end

% compute A matrix projected to POD basis
[Ur, Sr, Vr] = reducedsvd(Y0, r);
M = Y1*Vr*diag(1./diag(Sr));
tilA = Ur'*M;

% do eigendecomposition and project back eigenvectors
[tilw, lam, tilz] = eig(tilA);
lam = diag(lam);
w = M*tilw*diag(1./lam);
z = Ur*tilz;

% normalize the eigenvectors
for i=1:length(lam), z(:,i) = z(:,i)/(w(:,i)'*z(:,i)); end