function [lam, w, z] = subdmd(Y0, Y1, Y2, Y3, r)
% SUBDMD Compute subspace DMD.
%
% <Input>
%   Y0, Y1, Y2, Y3: N x M snapshot matrices, where N is the dimensionality
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
%   Takeishi et al., Subspace dynamic mode decomposition for stochastic
%     Koopman analysis, Phys. Rev. E 96:033310, 2017.
%     https://doi.org/10.1103/PhysRevE.96.033310

[n,m] = size(Y0);

% compute orthogonal projection of future onto past
Yp = [Y0; Y1]; Yf = [Y2; Y3];
[~, ~, Vp] = reducedsvd(Yp);
O = (Yf*Vp)*Vp';

% compute compact SVD of O and define Uq1 and Uq2
if nargin<5, r=rank(O); end
[Uq, ~, ~] = reducedsvd(O, r); r = min(n,r);
Uq1 = Uq(1:n,1:r);
Uq2 = Uq(n+1:end,1:r);

% compute projected A matrix
[U, S, V] = reducedsvd(Uq1);
M = Uq2*V*diag(1./diag(S));
tilA = U'*M;

% do eigendecomposition and project back eigenvectors
[tilw, lam, tilz] = eig(tilA);
lam = diag(lam);
w = M*tilw*diag(1./lam);
z = U*tilz;

% normalize the eigenvectors
for i=1:length(lam), z(:,i) = z(:,i)/(w(:,i)'*z(:,i)); end