% Second Harmonic Generation (SHG) Plane-Wave Analysis
% ------------------------------------------------

%% 1. Define Physical Constants and Parameters
clear; close all;
global eps0 c
c = 3e8;          % Speed of light (m/s)
eps0 = 8.85e-12;  % Vacuum permittivity (F/m)

% Material and field parameters
lambda1 = 1e-6;   % Fundamental wavelength (m)
n1 = 1.6;         % Refractive index at fundamental
deff = 1e-12;     % Effective nonlinear coefficient (m/V)
L = 1e-2;         % Medium length (m)
I1_0 = 1e13;      % Initial fundamental intensity (W/m²)

% Calculate nonlinear coupling coefficient eta
eta = 2*pi*c*deff/(n1*lambda1);

% Calculate initial field amplitude from intensity
A1_0 = sqrt(I1_0/(2*eps0*n1*c));

%% 2. Define Helper Functions

% Convert field amplitude to intensity
function I = get_intensity(A, n)
    global eps0 c
    I = 2*eps0*n*c*abs(A).^2;
end

% ODE system for SHG
function dA = shg_equations(z, A, params)
    eta = params.eta;
    dk = params.dk;
    
    dA = zeros(2,1);
    dA(1) = 1i*eta*A(2)*conj(A(1))*exp(-1i*dk*z);  % Fundamental
    dA(2) = 1i*eta*A(1)^2*exp(1i*dk*z);            % Second harmonic
end

%% 3. Phase-Matched Case (∆k = 0)
fprintf('Simulating phase-matched case...\n')

% Set up numerical solution
z = linspace(0, L, 1000);
params.eta = eta;
params.dk = 0;

% Solve coupled ODEs
[~, A] = ode45(@(z,A) shg_equations(z,A,params), z, [A1_0; 0]);

% Calculate intensities
I1_pm = get_intensity(A(:,1), n1);
I2_pm = get_intensity(A(:,2), n1);

% Plot results
figure(1)
plot(z*1000, I1_pm/I1_0, 'b-', 'LineWidth', 1.5)
hold on
plot(z*1000, I2_pm/I1_0, 'r--', 'LineWidth', 1.5)
xlabel('Distance (mm)')
ylabel('Normalized Intensity')
title('Phase-Matched SHG')
legend('Fundamental', 'Second Harmonic')
grid on

%% 4. Compare with Analytic Solution
fprintf('Comparing with analytic solution...\n')

% Calculate characteristic length
l = sqrt(2*n1^2*eps0*c^3/(2*pi*c/lambda1*deff*sqrt(I1_0)));

% Analytic solutions
A1_analytic = A1_0*sech(z/l);
A2_analytic = A1_0*tanh(z/l);

% Calculate intensities
I1_analytic = get_intensity(A1_analytic, n1);
I2_analytic = get_intensity(A2_analytic, n1);

% Plot comparison
figure(2)
plot(z*1000, I1_pm/I1_0, 'b-', z*1000, I1_analytic/I1_0, 'b--', ...
     z*1000, I2_pm/I1_0, 'r-', z*1000, I2_analytic/I1_0, 'r--', ...
     'LineWidth', 1.5)
xlabel('Distance (mm)')
ylabel('Normalized Intensity')
title('Numerical vs Analytic Solution')
legend('Numerical Fund.', 'Analytic Fund.', 'Numerical SH', 'Analytic SH', ...
       'Location', 'east')
grid on

%% 5. Phase-Mismatched Cases (Maker Fringes)
fprintf('Simulating phase-mismatched cases...\n')

% Reduced input intensity for better visibility
I1_0_reduced = 0.1e13;
A1_0_reduced = sqrt(I1_0_reduced/(2*eps0*n1*c));

% Different phase mismatches
dn_values = [1e-4, 3e-4];
figure(3)
hold on

for dn = dn_values
    % Calculate phase mismatch
    dk = 4*pi*dn/lambda1;
    params.dk = dk;
    
    % Solve ODEs
    [~, A] = ode45(@(z,A) shg_equations(z,A,params), z, [A1_0_reduced; 0]);
    I2 = get_intensity(A(:,2), n1);
    
    % Plot results
    plot(z*1000, I2/I1_0_reduced, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('\\Deltan = %.1e', dn))
    
    % Calculate and display coherence length
    Lcoh = pi/abs(dk);
    fprintf('Coherence length for Δn=%.1e: %.2f μm\n', dn, Lcoh*1e6)
end

xlabel('Distance (mm)')
ylabel('SH Conversion Efficiency')
title('Maker Fringes')
legend('show')
grid on

%% 6. Undepleted Pump Approximation
fprintf('\nSimulating undepleted pump approximation...\n')

% Function for undepleted pump integral
function A2 = undepleted_integral(z, A1_0, eta, dk)
    integrand = @(z_prime) A1_0^2 * exp(1i*dk*z_prime);
    A2 = 1i*eta*integral(integrand, 0, z);
end

% Calculate for specific phase mismatch
dn = 3e-4;
dk = 4*pi*dn/lambda1;

% Calculate SH field using integral
A2_undepleted = zeros(size(z));
for i = 1:length(z)
    A2_undepleted(i) = undepleted_integral(z(i), A1_0_reduced, eta, dk);
end
I2_undepleted = get_intensity(A2_undepleted, n1);

% Plot comparison
figure(4)
params.dk = dk;
[~, A] = ode45(@(z,A) shg_equations(z,A,params), z, [A1_0_reduced; 0]);
I2_numerical = get_intensity(A(:,2), n1);

plot(z*1000, I2_numerical/I1_0_reduced, 'b-', ...
     z*1000, I2_undepleted/I1_0_reduced, 'r--', ...
     'LineWidth', 1.5)
xlabel('Distance (mm)')
ylabel('SH Conversion Efficiency')
title('Numerical vs Undepleted Pump Approximation')
legend('Numerical', 'Undepleted Pump')
grid on