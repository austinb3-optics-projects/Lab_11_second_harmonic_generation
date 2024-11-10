% Function to calculate h(σ,ξ)
function h = calculate_h(sigma, xi)
    zeta = linspace(-xi, xi, 1000);
    dz = zeta(2) - zeta(1);
    integrand = exp(1i*sigma*zeta)./(1 + 1i*zeta);
    integral_result = sum(integrand)*dz;
    h = abs(integral_result)^2/(4.27*xi);
end

% Part (b) - h vs σ for fixed ξ
xi = 2.84;
sigma = linspace(-3, 3, 200);
h_sigma = zeros(size(sigma));
for i = 1:length(sigma)
    h_sigma(i) = calculate_h(sigma(i), xi);
end

figure(1);
plot(sigma, h_sigma);
xlabel('σ');
ylabel('h(σ,ξ)');
title(['Conversion Efficiency vs Phase Mismatch (ξ = ' num2str(xi) ')']);
grid on;

% Part (c) - h vs ξ for fixed σ
sigma_fixed = 0.56;
xi_range = linspace(1, 6, 200);
h_xi = zeros(size(xi_range));
for i = 1:length(xi_range)
    h_xi(i) = calculate_h(sigma_fixed, xi_range(i));
end

figure(2);
plot(xi_range, h_xi);
xlabel('ξ');
ylabel('h(σ,ξ)');
title(['Conversion Efficiency vs Focusing (σ = ' num2str(sigma_fixed) ')']);
grid on;

% Part (d) - 2D optimization
xi_2d = linspace(1, 4, 50);
sigma_2d = linspace(0, 2, 50);
[XI, SIGMA] = meshgrid(xi_2d, sigma_2d);
H = zeros(size(XI));

for i = 1:size(XI,1)
    for j = 1:size(XI,2)
        H(i,j) = calculate_h(SIGMA(i,j), XI(i,j));
    end
end

figure(3);
imagesc(xi_2d, sigma_2d, H);
colorbar;
xlabel('ξ');
ylabel('σ');
title('Conversion Efficiency h(σ,ξ)');
axis xy;

% Find optimal parameters
[max_h, idx] = max(H(:));
[opt_sigma_idx, opt_xi_idx] = ind2sub(size(H), idx);
opt_sigma = sigma_2d(opt_sigma_idx);
opt_xi = xi_2d(opt_xi_idx);

fprintf('Optimal parameters:\n');
fprintf('ξ_opt = %.2f\n', opt_xi);
fprintf('σ_opt = %.2f\n', opt_sigma);
fprintf('h_max = %.3f\n', max_h);

% Part (e) - Calculate physical parameters
L = 3e-3;  % Medium length (m)
lambda1 = 1e-6;  % Fundamental wavelength (m)
n1 = 1.6;  % Refractive index

% Calculate required spot size
b_opt = L/opt_xi;
w0_opt = sqrt(b_opt/(2*pi*n1/lambda1));

% Calculate required refractive index difference
dn_opt = 2*opt_sigma*lambda1/(pi*L);

fprintf('\nRequired physical parameters:\n');
fprintf('Spot size w₀ = %.2e m\n', w0_opt);
fprintf('Refractive index difference Δn = %.2e\n', dn_opt);

% Part (f) - Calculate conversion efficiency
deff = 15e-12;  % m/V
c = 3e8;  % m/s
eps0 = 8.85e-12;
P1 = 10;  % Input power (W)

K = 2.14*((2*pi*c/lambda1)^2)*deff^2*eps0*n1^2*n1/(c^3*pi);
ENL = K*L*(2*pi*n1/lambda1)*max_h;
P2 = ENL*P1^2;

fprintf('\nConversion parameters:\n');
fprintf('ENL = %.2e W⁻¹\n', ENL);
fprintf('Output power P₂ = %.2e W\n', P2);