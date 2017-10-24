function [Ur, Sr, Vr, contratio] = reducedsvd(X, r)
% REDUCEDSVD Compute reduced SVD of matrix X.

if nargin<2, r=rank(X); end

[Ur, Sr, Vr] = svd(X, 'econ');
diagSr = diag(Sr);
contratio = cumsum(diagSr)/sum(diagSr);
%r = min(r, sum(diagSr>eps));
Ur = Ur(:,1:r);
Sr = Sr(1:r, 1:r);
Vr = Vr(:,1:r);

end