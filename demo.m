% Demo script for subspace DMD. In this script, Koopman eigenvalues of a
% simple linear system (in two dimension) are compouted using subspace DMD
% and ordinary DMD (exact DMD).

rng(1234567890);

% settings
M = 1000;
noisestd_process = 0.1;
noisestd_observation = 0.1;
r = 0.9;
lam_true = [r*exp(1i*pi/180*(90)); r*exp(1i*pi/180*(-90))];

% generate data from a linear system
Y = zeros(length(lam_true), M);
Y(:,1) = ones(2,1);
for t=2:M
    Y(:,t) = diag(lam_true)*Y(:,t-1) + ...
        randn(length(lam_true),1)*noisestd_process;
end

% add observation noise
noise = randn(size(Y))*noisestd_observation;
[lam0, ~] = orddmd(Y(:,1:end-1), Y(:,2:end));
Y = Y + noise;

% subspace DMD
[lam_sub, w_sub, z_sub] = ...
    subdmd(Y(:,1:end-3), Y(:,2:end-2), Y(:,3:end-1), Y(:,4:end));
varphi_sub = z_sub'*Y; % values of Koopman eigenfunctions

% ordinary DMD
[lam_ord, w_ord, z_ord] = orddmd(Y(:,1:end-1), Y(:,2:end));
varphi_ord = z_ord'*Y;

% plot results
figure;
hold on;
plot(real(lam_true), imag(lam_true), 'o');
plot(real(lam_sub),  imag(lam_sub), 'o');
plot(real(lam_ord),  imag(lam_ord), 'o');
plot(sin(linspace(0,2*pi,100)), cos(linspace(0,2*pi,100)), 'k--');
hold off;
axis equal;
grid on;
xlabel('$\Re(\lambda)$', 'interpreter', 'latex');
ylabel('$\Im(\lambda)$', 'interpreter', 'latex');
legend({'truth','subspace DMD','ordinary DMD'});