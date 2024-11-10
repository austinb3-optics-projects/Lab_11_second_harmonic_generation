% Parameters
lambda1 = 1e-6;  % Fundamental wavelength (m)
n1 = 1.6;        % Refractive index
w0 = 10e-6;      % Beam waist (m)
k1 = 2*pi*n1/lambda1;  % Wave number
b = k1*w0^2;     % Confocal parameter

% Calculate z range (-3b to 3b)
z = linspace(-3*b, 3*b, 1000);

% Calculate w(z) and Φ(z)
w = w0*sqrt(1 + (2*z/b).^2);  % Spot size
phi = atan(2*z/b);            % Gouy phase

% Plot results
figure(1);
subplot(2,1,1);
plot(z/b, w/w0);
xlabel('z/b');
ylabel('w(z)/w₀');
title('Normalized Spot Size');
grid on;

subplot(2,1,2);
plot(z/b, phi/pi);
xlabel('z/b');
ylabel('Φ(z)/π');
title('Gouy Phase Shift');
grid on;